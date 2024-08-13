require(limSolve)
require(dplyr)
library(data.table)
library(readr)

# Function returning a list of window
# Input:
# size: the total size (all window cumul)
# winsize: the window size
# minsize: the minimum window size, i.e. the last window is smaller than winsize when size is not a multiple of winsize
#          => In this case, if the last window is smaller, fuse the two last windows
# Example:
# input: size = 45; winsize=10, minsize = 10
# output: list(1:10, 11:20, 21:30, 31:45)
get.win = function(size, winsize, minsize){
  
  if(winsize < size){
    
    nblocks = size/winsize
    fullblocks = floor(nblocks)
    
    ## split column numbers into 'nblocks' groups
    SPLIT <- split(1:(fullblocks*winsize), rep(1:fullblocks, each = winsize))
    if(nblocks>fullblocks){ 
      dblock = (nblocks-fullblocks)*winsize
      SPLIT[[length(SPLIT)+1]] = 1:dblock + max(fullblocks*winsize)
    }
    
    if(length(SPLIT[[length(SPLIT)]]) < minsize & length(SPLIT) > 1){
      SPLIT[[length(SPLIT)-1]] = c(SPLIT[[length(SPLIT)-1]],SPLIT[[length(SPLIT)]])
      SPLIT = SPLIT[1:(length(SPLIT)-1)]
    }
    
    
  }else{
    
    SPLIT = list(1:size)
  }
  
  return(SPLIT)
}



runlsei = function(gx,predictors){
  # gx = vector with single ril genotype with 0/0.5/1 for ref hom / het / alt hom
  # predictors in the founders genotype matrix for the same window
  # The function return a vector of length = ncol(predictors) with founder proportion
  # Example:
  # c(0,0,0,1) => 100% match with the 4th founder
  # c(0,0.2,0.8,0) => likely mean that there is a part of this genomic window corresponding to founder 1 and another to founder 2
  # c(0,0.5,0.5,0) can also mean that founder 2 and 3 are both an exact match (founders are identicals)
  
  # First, SNP with at least missing values
  predictNotMissing = apply(cbind(predictors, gx),1,function(x) sum(is.na(x))==0) 
  if(sum(predictNotMissing)<2){return(rep(NA, ncol(predictors)))}
  predictors = as.matrix(predictors[predictNotMissing,])
  gx = gx[predictNotMissing]
  
  #If a founder is heterozygous at some loci, consider it is the genotype of the ril at this loci
  # at SNP 10, ril is 0 and founder A is 0.5, change founder A SNP 10 in 0
  # It assumes that the founder is truely heterozygous (no genotyping error) 
  # At that whatever is the ril genotype at the herozygous snp should thus be a match
  # Or that we just don't know and still assume a match
  # You may want to question this choice because heterozygous allele is supposedly homozygous founder can be genotyping error
  # So you may just want change to founders[founders==0.5]=NA
  # If you do that a SNP which is heterozygous is one founder will be also considered as SNP with one missing value a be eliminated (above)
  # So you would loose information
  # Anyway, if the typical haplotype blocks contains a lot of different snp between founders it should not matters much
  # (If some founder are only different from a couple of SNP, it likely matters)
  heterozygous = which(predictors==0.5)
  if(length(heterozygous)>0){
    nsnp = nrow(predictors)
    rowhet=sapply(heterozygous, function(h){
      #wcol = floor(h/nsnp)+1
      wrow = (h/nsnp - floor(h/nsnp))*nsnp
      if(wrow == 0) wrow = nsnp#; wcol = wcol-1
      as.integer(round(wrow,digits = 1))
    })
    
    predictors[heterozygous]=gx[rowhet]
  }
  
  predictors = round(predictors)
  predictors[predictors>1]=1 
  
  # First look if there is a perfect match with a founder
  # If yes, the function return c(0,1,0,0) (perfect match with founder 2)
  # Alternatively, it would return something like  c(0.333, 0.333, 0.333, 0) if there is perfect match with multiple founders
  matchgeno = sapply(1:ncol(predictors), function(founder){ sum(predictors[,founder] == gx)/nrow(predictors)})
  matchgeno =  round(matchgeno,digits = 7) == 1
  if(sum(matchgeno)>0){out=ifelse(matchgeno,1/sum(matchgeno), 0); return(out)}
  
  #If not, use lsei to solves a least squares problem under both equality and inequality constraints
  # Basically, we get a vector of founder frequency that would explain well SNP frequencies in the ril
  d = ncol(predictors)
  A = predictors
  E = t(matrix(rep(1,d)))
  Ff = 1
  G = diag(rep(1,d))
  H = matrix(rep(0.000,d))
  Y = gx
  
  out = lsei(A=A,B=Y,E=E,F=Ff,G=G,H=H,verbose=F)
  
  if(out$IsError){out = rep(NA,d)}else{
    out = out$X
    
    # I don't think this part of the code is useful but it does noy harm. I will think about it.
    whichfounder = which(out>0)
    whichpolym = ispolymorphicSNP(subfounders=matrix(predictors[, whichfounder], ncol=length(whichfounder)))
    whichpolym = whichpolym[!(whichpolym %in% which(is.na(gx)))]
    matchgeno = sapply(whichfounder, function(founder){ sum(predictors[whichpolym,founder] == gx[whichpolym])})
    matchgeno = matchgeno==length(whichpolym) & length(whichpolym)>0
    
    if(sum(matchgeno)>0){y=rep(F,ncol(predictors)); y[whichfounder]=matchgeno; out=ifelse(y,1/sum(y), 0)}
  }
  
  return(out)
  
}



pushleft = function(b1,b2,ril, founders, nsnp, rule = 1){
  # This function start from a window between b1 & b2
  # and increase the window range on the left (=decrease b1)
  # to check what how it affects the output of runlsei
  # And to stop once optimize the chosen parameter
  # rule 1 = maximize the max founder weight
  # rule 2 = minimize the number of match founder
  # Basically, this function is useful if we know that window b1=>b2 does match with founder X
  # and that we want to know the left boudary of this haplotype block
  # If at b1 - 10 the output of runlsei still indicate a perfect match with founder X, 
  # it indicates that the founder X haplotype block left boundaries can be pushed by -10
  # if the founder weights output by lsei decrease, it mean that we cannot push the boundary
  # i.e. the breakpoint is before
  
  b1push = b1
  minsize = nsnp + 1
  ni=0
  
  outlsei = runlsei(gx=ril[b1:b2],predictors = founders[b1:b2,])
  maxw = round(max(outlsei), digits = 3)
  nmatch = vtar = length(which(outlsei>0))
  
  if(rule==1){
    vtar = maxw # target value (max weight)
    # condition = continue to "push" the border while we increase the weight
    condition = (maxw >= vtar & b1push > minsize)
  } 
  
  if(rule==2){
    vtar = nmatch # target value (the number of founder with weight > 0 )
    # condition = continue to "push" the border while we reduce or keep the same number of founder with a match 
    condition = (nmatch <= vtar & b1push > minsize) 
  } 
  
  if(is.na(vtar)){stop()}
  
  while(condition==T){
    ni = ni+1
    b1push = b1push-nsnp
    win = b1push:b2
    predictors = founders[win,]
    gx = ril[win]
    
    outlsei = try(runlsei(gx=gx,predictors = predictors),silent = T)
    maxw = round(max(outlsei), digits = 3)
    nmatch = length(which(outlsei>0))
    
    if(rule == 1){
      if(is.na(maxw)) maxw = 0 # some trick when the haplotype weight search fail
      # condition = continue to "push" the border while we increase the weight
      condition = (maxw >= vtar & b1push > minsize)
      vtar = maxw
    }
    
    if(rule==2){
      if(nmatch==0) nmatch = 1e6 # some trick when the haplotype weight search fail
      # condition = continue to "push" the border while we minimize the number of founder with a match
      condition = (nmatch <= vtar & b1push > minsize)
      vtar = nmatch 
    } 
    
    
  }
  
  b1push = b1push+ifelse(ni>0, nsnp, 0)
  return(b1push)
}



pushright = function(b1,b2,ril, founders, nsnp, rule = 1){
  #rule 1 = maximaize the max weight
  #rule 2 = minimize the number of match founder
  # see pushleft function for details
  
  b2push = b2
  maxsize = length(ril)-nsnp
  ni=0
  
  outlsei = runlsei(gx=ril[b1:b2],predictors = founders[b1:b2,])
  maxw = round(max(outlsei), digits = 3)
  nmatch = vtar = length(which(outlsei>0))
  
  if(rule==1){
    vtar = maxw
    condition = (maxw >= vtar & b2push < maxsize)
  } 
  
  if(rule==2){
    vtar = nmatch 
    condition = (nmatch <= vtar & b2push < maxsize)
  } 
  
  if(is.na(vtar)){stop()}
  
  while(condition==T){
    ni = ni+1
    b2push = b2push+nsnp
    win = b1:b2push
    predictors = founders[win,]
    gx = ril[win]
    
    outlsei = try(runlsei(gx=gx,predictors = predictors),silent = T)
    maxw = round(max(outlsei), digits = 3)
    nmatch = length(which(outlsei>0))
    
    if(rule == 1){
      if(is.na(maxw)) maxw = 0
      condition = (maxw >= vtar & b2push < maxsize)
      vtar = maxw
    }
    
    if(rule==2){
      if(nmatch==0) nmatch = 1e6
      condition = (nmatch <= vtar & b2push < maxsize)
      vtar = nmatch 
    } 
    
    
  }
  
  b2push = b2push-ifelse(ni>0, nsnp, 0)
  return(b2push)
}


# pushright = function(b1,b2,ril, founders, nsnp){
#   
#   maxw = vtar = round(max(runlsei(gx=ril[b1:b2],predictors = founders[b1:b2,])), digits = 3) 
#   if(is.na(vtar)){stop()}
#   b2push = b2
#   
#   maxsize = length(ril)-nsnp
#   ni = 0
#   
#   while(maxw >= vtar & b2push < maxsize){
#     ni = ni+1
#     b2push = b2push+nsnp
#     win = b1:b2push
#     predictors = founders[win,]
#     gx = ril[win]
#     
#     maxw = try(max(runlsei(gx=gx,predictors = predictors)))
#     if(is.na(maxw)){maxw=0}
#     maxw = round(maxw, digits = 3)
#   }
#   
#   b2push = b2push-ifelse(ni>0, nsnp, 0)
#   return(b2push)
# }


ispolymorphicSNP = function(subfounders){
  #return where are the polymorphic snps among founders (generally subset of all founders)
  if(ncol(subfounders)>1){
    mx = apply(subfounders, 1, mean)
    dx = apply(subfounders,2,function(x){abs(x-mx)})
    whichpolym = which(apply(dx,1,sum)>0)
  }else{
    whichpolym=which(F)
  }
  
  return(whichpolym)
}




find.border = function(b1,b2,ril, founders,rule=1){
  # Starting from a initial window between b1 & b2 matching a founder
  # the function expand (push) the window to find the boundaries/breakpoints of the founder haplotype block
  # It use pushleft and pushright function (more comments in pushleft function)
  b1push=pushleft(b1=b1,b2=b2,ril=ril, founders=founders, nsnp=50, rule=rule)
  b1push=pushleft(b1=b1push,b2=b2,ril=ril, founders=founders, nsnp=10, rule=rule)
  b1=pushleft(b1=b1push,b2=b2,ril=ril, founders=founders, nsnp=1, rule=rule)
  
  
  #Search for right border
  b2push = pushright(b1=b1,b2=b2,ril=ril, founders=founders, nsnp=50, rule=rule)
  b2push = pushright(b1=b1,b2=b2push,ril=ril, founders=founders, nsnp=10, rule=rule)
  b2 = pushright(b1=b1,b2=b2push,ril=ril, founders=founders, nsnp=1, rule=rule)
  
  # Return the founder haplotype
  win = b1:b2
  predictors = founders[win,]
  gx = ril[win]
  
  wfounder = runlsei(gx=gx,predictors = predictors)
  
  return(cbind(data.frame(pos1 = b1, pos2=b2), matrix(wfounder, nrow=1)))
  
}

supressEncased = function(HAPLO){
  
  # Search for inferred haplotype encased in another
  # i.e. hap1 go from 1 to 150 & hap2 go from 10 to 140
  # and choose the larger haplo
  
  HAPLO = as.data.frame(HAPLO)
  HAPLO = HAPLO[order(HAPLO$pos1,HAPLO$pos2),]
  
  h = 1
  while(h < nrow(HAPLO)){
    
    pos1=HAPLO[h,1];pos2=HAPLO[h,2]
    
    isEncased = HAPLO[,1] <= pos1 & HAPLO[,2] >= pos2
    if(sum(isEncased)>1){HAPLO=HAPLO[-h,]}else{h=h+1}
  }
  
  return(HAPLO)
}


splithaplo = function(b1,b2,ril, founders){
  
  # when the output of runlsei indicate that between b1 & b2 the haplotype is a combination of several founders
  # This particular function can be used to try to split the window to find the breakpoint between  the different founder haplotypes
  # i.e. within win to 1 to 100 there is a combination of founder A and B
  # the function can return that 1 to 70 is founder B and 71 to 100 is founder A
  
  # First we want to find which are the polymorphic snps between the possible founders
  win = b1:b2
  gx = ril[win]
  predictors = as.matrix(founders[win,])
  
  founderweights = runlsei(gx,predictors)
  whichfounder =  which(founderweights>0)
  
  heterozygous = which(predictors==0.5)
  if(length(heterozygous)>0){
    nsnp = nrow(predictors)
    rowhet=sapply(heterozygous, function(h){
      #wcol = floor(h/nsnp)+1
      wrow = (h/nsnp - floor(h/nsnp))*nsnp
      if(wrow == 0) wrow = nsnp#; wcol = wcol-1
      as.integer(round(wrow,digits = 1))
    })
    
    predictors[heterozygous]=gx[rowhet]
  }
  
  whichpolym = ispolymorphicSNP(subfounders=predictors[, whichfounder])
  whichpolym = whichpolym[!(whichpolym %in% which(is.na(gx)))]
  
  #View(cbind(gx, predictors[,whichfounder]))
  
  # Then, at the polymorphic snps, we look at the that corresponding founder
  # search for indication of breakpoints between the founders
  if(length(whichpolym)>0){
    
    breaks = c(0,whichpolym,length(win)+1)#+b1-1
    breaks = t(sapply(2:(length(breaks)-1), function(i){c(breaks[i-1]+1, breaks[i+1]-1)}))
    breaks = matrix(breaks, ncol=2)
    matchgeno = sapply(whichpolym, function(wsnp){ which(predictors[wsnp,whichfounder] == gx[wsnp])})
    matchgeno = matrix(matchgeno, nrow = length(whichpolym))
    matchgeno = t(sapply(whichpolym, function(wsnp){
      out = rep(NA, length(whichfounder))
      wx = which(predictors[wsnp,whichfounder] == gx[wsnp])
      out[wx]=1
      out[is.na(out)]=0
      out}))
    
    if(nrow(matchgeno)>1){
      matchgeno=apply(matchgeno, 2, function(x){
        wbreak = diff(c(0,which(x==0), length(x)+1))-1
        wbreak = wbreak[wbreak>0]
        wbreak = unlist(sapply(wbreak, function(wb){rep(wb,wb)}))
        x[x==1]=wbreak
        x
      })
    }
    
    
    matchgeno2 = apply(matchgeno, 1, which.max)
    
    tofuse = which(diff(matchgeno2)==0)
    if(length(tofuse)>0){
      for(tf in tofuse){breaks[tf+1,1] = breaks[tf,1]}
      breaks = breaks[-tofuse,]
    }
    
    breaks = matrix(breaks, ncol=2)
    #1:nrow(breaks)
    splitted = do.call(rbind,lapply(1:nrow(breaks), function(i){
      #print(i)
      start = breaks[i,1]
      end = breaks[i,2]
      outf = try(runlsei(gx=gx[start:end],predictors=predictors[start:end,]),silent = T)
      if(class(outf)=="try-error") outf = matrix(rep(NA, ncol(predictors)), nrow=1)
      
      x = cbind(data.frame(pos1 = start+b1-1, pos2=end+b1-1), matrix(outf, nrow=1))
    }))
    
    return(splitted)
    
  }else{
    return(cbind(data.frame(pos1 = b1, pos2=b2), matrix(founderweights, nrow=1)))
  }
  
  
  
}





haplosearch = function(ril, founders, snps){
  ############################
  #### Main function #########
  ############################
  
  # 1) First we try to find regions with good match with some founders
  # We look within many different windows with varying size for correspondence with founders
  # This is done with runlsei out with return a vector of weight of lenght = ncol(founders)
  # 0,1,0,0 Perfect match with founder 2
  # several numbers > 0 & < 1, several founders possible or combination of founders
  # see runlsei for details
  winsizes = c(300, 200, 150, 125, 100, 75, 50)
  # For each window size:
  lsout = lapply(winsizes, function(wsize){
    
    # Get a list windows along the chromosome
    windows = get.win(nrow(founders), winsize=wsize, minsize=50)
    
    # For each window, test for match with founder haplotype using the runlsei function
    out = do.call(cbind, lapply(windows, function(win){
      
      predictors = founders[win,] # predictors = founders in the window
      gx = as.numeric(ril[win]) # gx = ril genotype in the window
      
      # get founder weights for this window
      p = try(runlsei(gx=gx,predictors = predictors), silent = T)
      if(class(p)=="try error"){p = rep(NA,ncol(predictors))}
      
      p=as.numeric(p)
     
      # Repeat the weights over window size (at the end we have a matrix with the weight for each snp)
      p=matrix(rep(p,  length(win)), ncol=length(win), nrow=ncol(predictors), byrow = F)
      p
      
    }))
    
    out
  })
  
  names(lsout)=winsizes

  
  # lsout = matrix with row correponding to different window size
  # columns correspoinding to snps
  # and value to the max founder weights 
  # Select the max value
  # = select the window size a each position that maximize the founder weights (ideally a perfect match to one founder = weight of 1)
  maxf = do.call(rbind, lapply(lsout, function(x){apply(x,2,max)})) # matrix with max founder weight value, row = different winsizes, column = different snps
  maxf = unlist(apply(maxf, 2, function(x){ if(sum(is.na(x))<nrow(maxf)){which.max(x)}else{NA}})) # vector indicating which window size maximize the max founder weight for each snp
  
  # Now, use maxf to extract the founder weights, chosing for each snp the window size maximizing the founder weight
  # The weights inform us about the match with the different founders 
  # haploweights is a matrix with row = foudners; col = snps;
  # and values are the founder' weights (chosen among the different window to maxize the max weight)
  haploweights = do.call(cbind, lapply(1:length(maxf), function(i){
    if(!is.na(maxf[i])){as.numeric(lsout[[ maxf[i] ]][,i])}else{rep(NA, ncol(founders))}
  }))
  
  haploweights = round(haploweights, digits = 3)
  # haploweights inform us about regions that have a good match with some founder
  # i.e., haploweight[i:j] = 1, means that between i:j there is a perfect match with a founder
  # We want to find the limit of the corresponding founder haplotype block (here I'll call these limits breakpoints)

  
  
  maxweight = apply(haploweights, 2, max)
  
  # Now that we have tested for different windows along the chromosome
  # and kept the ones with the maximum max weights (i.e., with the highest match with a founder)
  # => stored in haploweights (weights for all founders) and maxweight (only the max founder weights)
  # We know regions that have a good match with a founder (ideally, a perfect one)
  # We want to find the limit of these regions, i.e. the breakpoints of an haplotypes (note: these limits will be refined later)
  # We are going to extract the breakpoints of haplotype starting with the ones with the highest maximum weights (best match with a founder)
  # Then refine the haplotype limit and optimize the match with a single founders. 
  # Results will be stored in HAPLO
  ### Note: we start with the highest weight window because we will update the window by refining the limits
  ### and some most weight window (region that are not specific to a founder) will be incorporated in other haplotypes when expending limits, so we won't have to deal with them
  uniqmaxweight = sort(unique(maxweight[maxweight>0]),decreasing = T) # All the unique weights
  HAPLO = NULL
  
  while(max(uniqmaxweight)>0){
    
    #Target value: the maximum value of uniqmaxweight vector
    vtar = max(uniqmaxweight,na.rm=T)
    #print(vtar)
    
    # Which snp have the vtar value of maxweight
    wvtar = which(maxweight==vtar)
    
    # wvtar can contain multiple primer of haplo
    
    # What I call breakpoints are the recombinations points between haplotypes
    # Even though at this stage their position is unprecise (not have been refined)
    # We want to infer these unprecise breakpoints from 
    # Here, breakpoints is a two column matrix with the wvtar vector
    # i.e., where there is a gap between snps. 
    # 1 2 3 4 10 11 12 14 => breakpoints = 1,4,10,14
    
    breakpoints = sort(c(wvtar[1], # first snp position of wvtar (1 in our example above)
                         wvtar[which(diff(wvtar)>1)], # The snp before the gap(s) (4 in out example aove)
                         wvtar[which(diff(wvtar)>1)+1], #the snp after the gap(s) (10 in out example aove)
                         wvtar[length(wvtar)])) # The last snp (14 in out example aove)
    
    # transform a matrix where the each row = c(breakpoint 1, breakpoints 2) = the limits of one haplotype
    breakpoints = matrix(breakpoints, ncol=2, byrow=T) 
    # # I guess I put this line to resolve a bug some unusual case where the brekpoints are not in the good order
    breakpoints = matrix(breakpoints[breakpoints[,2]-breakpoints[,1] > 0,], ncol=2) 
    
    if(nrow(breakpoints)>0){
      
      # For now, this approach does not ensure that a line of breakpoints correspond to a single haplotype
      ### i.e., 1:1000 have a perfect match with some founder (weight of 1) => breakpoints[i,] = c(1,1000)
      ### but 1:600 can corresponds to founder 1 and 601:1000 to founder 2
      # => The next block of code checks this possibility and add breakpoints when needed (i.e. c(1,1000) => c(1,600) & c(601,1000))
      breakpoints = do.call(rbind, lapply(1:nrow(breakpoints), function(bi){
        b1 = breakpoints[bi,1] #breakpoints 1 ( left limit of the haplotype corresponding to breakpoints[bi,])
        b2 = breakpoints[bi,2] #breakpoints 2 ( right limit of the haplotype corresponding to breakpoints[bi,])
        
        # Look between the breakpoints (b1:b2) which founder have the maximum weights at each snp => store in wfounder
        # If different founders in the target region, n = length(unique(wfounder)) is greater than 1
        wfounder = apply(haploweights[,b1:b2], 2, which.max)
        n = length(unique(wfounder))
        
        # if more than one founder add split  breakpoints[bi,] in two adding the new breakpoints
        if(n>1){
          newbreaks = b1 + which(diff(wfounder)!=0)
          #newbreaks = b1+cumsum(unlist(lapply(split(wfounder,wfounder), length)))
          newbreaks = c(b1,newbreaks,b2+1)
          newbreaks = cbind(newbreaks[1:(length(newbreaks)-1)], newbreaks[2:length(newbreaks)]-1)
          return(newbreaks)
        }else{
          return(breakpoints[bi,])
        }
        
      }))
      
      
      # Now we have the pre haplotype blocks limits stored in brekpoints
      ### i.e., we founder that snp100:snp1000 had a perfect match with founder 1
      ###       but because we did not tested all the possible window ,
      ###       it is very likely that the limits are not 100 & 1000 but something like 75 and 1093
      ###       Also for regions that don't have a prefect match (multiple weights > 0 & < 1),
      ###       it could mean that the region contain a singl haplotype but e is not specific to a single founder
      ###       or the region contains multiple haplotype matching with different founder
      ###       If it is the case, we're going to try split the haplotype to see if we can find multiple haplo with a better match
      haplo=do.call(rbind,lapply(1:nrow(breakpoints), function(bi){
        
        #For each (pre) haplotype blocks
        # Extracts brekpoints (haplotype block limits)
        b1 = breakpoints[bi,1]
        b2 = breakpoints[bi,2]
        
        # Use the function find.border to find the real haplotypes limits/brekpoints
        # i.e, expand the limit to farthest possible while optimizing the match with a founder
        # See find.border function for more details
        x = try(find.border(b1=b1,b2=b2,ril=ril, founders=founders),silent = T)
        if(class(x)=="try-error") x = matrix(c(b1,b2,rep(NA, ncol(founders))), nrow=1)
        
        # x is like this:
        # pos1 pos2 1 2 3 4
        #   3   819 0 0 0 1
        # => pos1 and pos2 are the haplotype limits/breakpoints,
        #    and 1,2,3,4,... are the founder weights indicating the different weights 
        
        # which founders have weight > 0
        wfounder = which(x[,3:ncol(x)]>0)
        
        # length(wfounder > 1): no prefect match with a single founder (multiple weights > 0 & < 1),
        ### It could mean that the region contain a single haplotype but is not specific to a single founder
        ### or the region contains multiple haplotype matching with different founder
        # =>  Using splithaplo, we're going to try split the haplotype to see if we can find multiple haplo with a better match
        # The splithaplo function try to split the haplo is a 'wise' way focusing of the existing polymorphism between the possible founders
        
        if(length(wfounder)>1){
          pos1 = x[,1]
          pos2 = x[,2]
          
          x=splithaplo(b1=pos1,b2=pos2,ril=ril, founders=founders)
        }
        colnames(x) = c("pos1", "pos2", 1:ncol(founders))
        return(x)
        
        
      }))
      
      
      colnames(haplo)=c("pos1", "pos2", 1:ncol(founders))
      
      # haplo look like this
      #pos1 pos2 1 2 3 4
      #  3  819  0 0 0 1
      # 760 2001 1 0 0 0
      # => each row is one haplotype block, pos1/pos2 its limits (recombination breakpoints) and the match'weight with the different founders
      ###  Note: in this example the two haplotype block overlaps. It means that between snp 760 and 819 founder 1 and 4 have the same sequence
      ###        So we don't know exactly where the recombination breakpoints happened between those two positions 
      
      # Here, we "update" haploweight matrix (the one which contains the initial weights from the multiple window we tested)
      # Puting a -1 where are the haplotypes blocks we just infered
      # This is just a trick so we won't check these positions a second time
      # Ex:
      # initial haplomatrix
      # 1 1 1 1 1 1 1 1 1 0.5 0.5 0.5 0.5 0.5
      # 0 0 0 0 0 0 0 0 0 0.5 0.5 0.5 0.5 0.5
      # 0 0 0 0 0 0 0 0 0  0   0   0   0  0
      # => by extending the limit of the haplotype with weight == 1 for founder 1, we found that the whole region corresponds to founder 1
      # So we don't want to lose time trying to find the limit of the initial block with weights == 0.5 for founder 1 and 2
      # We're just going to replace all the value by -1 so we know that we don't need to investiate them again
      # Note: the founder weight obtained are NA, we replace the value by NA rather than -1
      for(bi in 1:nrow(haplo)){ # for each inferred founder haplo block
        # the limits/brekpoints
        b1 = haplo[bi,1] 
        b2 = haplo[bi,2]
        # replace haploweights[,b1:b2] by -1 (or NA)
        v = ifelse(is.na(haplo[bi,3]), NA, -1)
        fx = matrix(v, nrow = ncol(haplo)-2, ncol = length(b1:b2))
        haploweights[,b1:b2] = fx
      }
      
      # Last, bind haplo (haplotype found for this round of search) to HAPLO (all the haplo we find)
      HAPLO = rbind(HAPLO, haplo)
      
    }else{ 
      haploweights[,wvtar] = NA 
    }
    
    # Update unique maxweight
    # As we just changed the weight value to -1 (or NA) for the initial highest weight we investigated
    # The second initial highest weight is now ranked first
    maxweight = apply(haploweights, 2, max)
    uniqmaxweight = sort(unique(maxweight[maxweight>0]),decreasing = T)
    
  }
  
  # Get rid of NA
  HAPLO = HAPLO[!is.na(HAPLO[,3]),]
  
  # HAPLO look like this:
  # pos1 pos2 1 2 3 4
  #   3  819 0 0 0 1
  #  760 2001 1 0 0 0
  #   1   819 0 0 0 1
  # => at this stage: some haplotype overlap, some need to be fused
  
  # But first, we redoing a second try to expand the border 
  # It is likely needed because of the haplotype we split, but not sure
  # Worst case scenario, this double check is not useful and we just slow down a bit the code
  HAPLO  = do.call(rbind, lapply(1:nrow(HAPLO), function(i){
    #print(i)
    pos1 = unlist(HAPLO[i,1])
    pos2 = unlist(HAPLO[i,2])
    
    x = find.border(b1=pos1,b2=pos2,ril=ril, founders=founders, rule=1)
    x = round(unlist(x), digit = 5)
    x = matrix(x, nrow=1)
    if(length(which(x[,3:ncol(x)]>0))>1){x=splithaplo(b1=x[,1],b2=x[,2],ril, founders)}
    
    colnames(x) = colnames(HAPLO)
    
    x
  }))
  
  HAPLO = HAPLO[!is.na(HAPLO[,3]),]
  
  
  # As seen above, in some haplotype blocks stored in HAPLO overlap (or are the same, found two times from different initial window)
  # if one haplo is between 1 and 1000 and the second between 200 and 900
  # Just consider the first one
  HAPLO = supressEncased(HAPLO)
  # Now HAPLO look like this:
  # pos1 pos2 1 2 3 4
  #  1  819   0 0 0 1
  # 760 2001  1 0 0 0
  
  
  
  # Next fuse overlaping windows with common founders
  # haplotype 1 is between 200 and 1000 + match to founder 1; haplotype 2 between 600 and 1200 + match to founder 1 and 2 (not specific)
  # Fuse them => New haplotype between 200 and 1200, match to founder 1
  # I think these case happen from the split haplotype (or because find.border in not perfect and sometimes does not push the border enough?)
  # I would need to dive in to that again to clarify extactly when these occurs
  # Or if this part of the code the remains of the draft code
  # But for now, that's do the job
  
  HAPLO = HAPLO[order(HAPLO$pos1,HAPLO$pos2),]
  
  i = 1
  while(i < nrow(HAPLO)){
    j=i+1
    #pos1=HAPLO[i,1];pos2=HAPLO[i,2]
    overlapwithnext = HAPLO[j,1] <= HAPLO[i,2]
    
    x1 = HAPLO[i,3:ncol(HAPLO)]
    x2 = HAPLO[j,3:ncol(HAPLO)]
    x1[x1>0] = round(1, digits=0)
    x2[x2>0] = round(1, digits=0)
    
    commonfounder = which((x1 + x2) > 1)
    
    if(length(commonfounder)>0 & overlapwithnext){
      b1 = HAPLO[i,1]
      b2 =HAPLO[j,2]
      x = round(unlist(runlsei(ril[b1:b2],founders[b1:b2,])), digits = 2)
      
      if(!is.na(x[1])){
        x = c(b1,b2,x)
        HAPLO[i,] = x
        HAPLO = HAPLO[-j,]
      }
    }
    
    i = i+1
    
    #}else{i = i+1}
  }
  
  HAPLO = HAPLO[!is.na(HAPLO[,3]),]
  
  HAPLO = HAPLO[order(HAPLO$pos1,HAPLO$pos2),]
  
  colnames(HAPLO) = c("whichsnp1", "whichsnp2", colnames(founders))
  
  return(HAPLO)
  
}





haploDfToMatrix = function(haplotype, snps){
  
  haplomatrix = lapply(split(haplotype, haplotype$rilname), function(x){
    print(x$rilname[1])
    #x=subset(haplotype,rilname == 101)
    xx = do.call(rbind,lapply(1:nrow(x), function(i){
      #print(i)
      y=x[i,] 
      w1 = which.min(abs(y$pos1-snps$pos))
      w2 = which.min(abs(y$pos2-snps$pos))
      founder = which(y[,3:(ncol(y)-1)]>0)
      if(length(founder)>0){
        return(cbind(w1:w2, rep(founder, w2-w1+1)) )}else{
          return(NULL)
        }
    }))
    
    mm = match(1:max(xx[,1]), xx[,1])
    xx = xx[mm,2]
    xx = matrix(xx, ncol=1)
    colnames(xx) = x$rilname[1]
    xx
    #unlist(xx)
  })
  
  nsnp = max(unlist(lapply(haplomatrix, nrow)))
  
  haplomatrix = do.call(cbind, lapply(haplomatrix,function(x){ x[match(1:nsnp,1:nrow(x))] }))
  
  haplomatrix
}




KeepMostFrequentHaplo = function(haplotype, rils){
  hapfreq = do.call(rbind, lapply(1:nrow(rils), function(i){
    if(i %% 1000 == 0){print(i)}
    pos = snps$pos[i]
    hap = haplotype[haplotype[,1] <= pos & haplotype[,2] > pos,]
    apply(hap[,3:(ncol(hap)-1)], 2, sum)/nrow(hap)
    
  }))
  
  #haplotype haplotype2
  haplotype = do.call(rbind, lapply(1:nrow(haplotype), function(i){
    if(i %% 1000 == 0){print(i)}
    x = haplotype[i,]
    wfounder = x[,3:(ncol(x)-1)]>0
    if(sum(wfounder)>1){
      pos1 = x[,1]
      pos2 = x[,2]
      pos1 = which.min(abs(snps$pos-pos1))
      pos2 = which.min(abs(snps$pos-pos2))
      
      freq = hapfreq[pos1:pos2,]#[wfounder]
      ismax = apply(freq, 2, mean)[wfounder]
      ismax = which.max(ismax)
      wfounder = which(wfounder)[ismax]
      x[,3:(ncol(x)-1)] = 0
      x[,wfounder+2] = 1
    }
    
    return(x)
  }))
  
}



hapstat = function(haplotype, winsize=10000){
  haplostat = do.call(rbind, lapply(seq(1, max(haplotype[,2]), winsize), function(pos){
    
    keep = haplotype[,1] <= pos & haplotype[,2] > pos
    x=haplotype[keep,]
    haplength = x[,2] - x[,1]
    nhap = sum(apply(x[,3:(ncol(x)-1)],2, sum) > 0)
    
    data.frame(pos = pos,
               mean.size = mean(haplength, na.rm=T),
               min.size = min(haplength, na.rm=T),
               max.size = max(haplength, na.rm=T),
               nhap = nhap)
    
  }))
  
  haplostat
}











FuseHapSeparatedByOneSNP = function(HAPLO, posinfo = c("whichsnp", "pos", "cM"), foundercolnames=c("FA.g1", "FA.g2", "FM.g1", "FM.g2")){
  
  HAPLO = HAPLO[order(HAPLO$whichsnp1,HAPLO$whichsnp2),]
  
  fused = do.call(rbind, lapply(1:nrow(HAPLO), function(i){
    #print(i)
    whichfounder = which(HAPLO[i,foundercolnames]>0)
    founder = matrix(unlist(round(HAPLO[,whichfounder+2])), nrow = nrow(HAPLO))
    founder[founder>0.0] = 1
    
    #Look if there is another haplotype block with corresponding to the same founder
    nextcommon =  sort(unlist( apply(founder, 2, function(x){which(x>0)}) ))
    nextcommon = min(nextcommon[nextcommon >i])
    if(nextcommon==Inf){return(NULL)}
    j = nextcommon 
    
    #Look if there is an overlap between haplotype i & j
    overlapwithnext = HAPLO[j,1] <= HAPLO[i,2]
    
    #gap size between i and j
    gapsize = HAPLO[j,"whichsnp1"] - HAPLO[i,"whichsnp2"] - 1
    
    if(!overlapwithnext & gapsize == 1){
      
      commonfounder = HAPLO[i,foundercolnames] > 0 & HAPLO[j,foundercolnames] > 0
      commonfounder = as.numeric(commonfounder)/sum(commonfounder)
      new = HAPLO[i,]
      new[,foundercolnames] = commonfounder
      new$whichsnp2 = HAPLO[j,]$whichsnp2
      if("pos" %in% posinfo) new$pos2 = HAPLO[j,]$pos2
      if("cM" %in% posinfo) new$cM2 = HAPLO[j,]$cM2
      return(new)
      
    }else{return(NULL)}
  }))
  
  if(is.null(fused)){fused = HAPLO[which(F),]}
  ### Fuse the fused haplotype that overlap
  i = 1
  while(i < nrow(fused)){
    j=i+1
    #pos1=HAPLO[i,1];pos2=HAPLO[i,2]
    overlapwithnext = fused[j,"whichsnp1"] <= fused[i,"whichsnp2"]
    commonfounder = fused[i,foundercolnames]>0 & fused[i,foundercolnames]>0
    
    if(sum(commonfounder)>0 & overlapwithnext){
      new = fused[i,]
      new[,foundercolnames] = commonfounder
      new$whichsnp2 = fused[j,]$whichsnp2
      if("pos" %in% posinfo) new$pos2 = fused[j,]$pos2
      if("cM" %in% posinfo) new$cM2 = fused[j,]$cM2
      
      fused[i,]=new
      fused=fused[-j,]
    }else{i = i+1}
  }
  
  
  
  
  HAPLO = as.data.frame(rbind(HAPLO, fused))
  
  HAPLO = HAPLO[order(HAPLO$whichsnp1,HAPLO$whichsnp2),]
  
  h = 1
  while(h < nrow(HAPLO)){
    
    snp1=HAPLO[h,"whichsnp1"];snp2=HAPLO[h,"whichsnp2"]
    
    isEncased = HAPLO[,1] <= snp1 & HAPLO[,2] >= snp2
    if(sum(isEncased)>1){HAPLO=HAPLO[-h,]}else{h=h+1}
  }
  
  HAPLO
}



putBreakAtMidGeneticDistance=function(HAPLO,info){
  HAPLO = HAPLO[order(HAPLO$whichsnp1,HAPLO$whichsnp2),]
  
  for(i in 1:(nrow(HAPLO)-1)){
    j = i + 1
    
    if(HAPLO[i,"whichsnp2"]>HAPLO[j,"whichsnp1"]){
      cMi = HAPLO[i,"cM2"]
      cMj = HAPLO[j,"cM1"]
      
      newcM = (cMi+cMj)/2
      newpos = approx(info$cM,info$pos, xout = newcM)$y
      newsnp = which.min(abs(info$cM-newcM))
      
      HAPLO[i,"cM2"] = newcM
      HAPLO[j,"cM1"] = newcM
      
      HAPLO[i,"pos2"] = newpos
      HAPLO[j,"pos1"] = newpos
      
      HAPLO[i,"whichsnp2"] = newsnp-1
      HAPLO[j,"whichsnp1"] = newsnp
    }
  }
  
  #HAPLO = HAPLO[,-which(colnames(HAPLO) %in% c("whichsnp1", "whichsnp2"))]
  HAPLO
}



divideOverlapHapBlocksByFounderID = function(haplotype, foundercolnames=c("FA.g1", "FA.g2", "FM.g1", "FM.g2"), info){
  
  haplotype=haplotype[order(haplotype$whichsnp1, haplotype$whichsnp2),]
  
  snpmin = min(haplotype$whichsnp1)
  snpmax = max(haplotype$whichsnp2)
  
  rilname=haplotype$rilname[1]
  
  snphapmatrix=do.call(rbind, lapply(foundercolnames, function(thisfounder){
    
    whichsnphap = 1:snpmax
    snphap = rep(0, length(whichsnphap))
    thisfounderhaplotypes=haplotype[haplotype[,thisfounder]>0,]
    if(nrow(thisfounderhaplotypes)>0){
      for(i in 1:nrow(thisfounderhaplotypes)){
        wsnp1 = thisfounderhaplotypes[i,]$whichsnp1
        wsnp2 = thisfounderhaplotypes[i,]$whichsnp2
        snphap[wsnp1:wsnp2]=1
      }
    }
    
    snphap
  }))
  
  
  snphapid=apply(snphapmatrix, 2, function(x){paste(foundercolnames[x>0], collapse=";")})
  
  breaks = which(diff(as.numeric(factor(snphapid)))!=0)
  breaks=breaks[breaks>snpmin]
  breaks = c(snpmin, breaks, snpmax)
  
  
  haplotype=do.call(rbind,lapply(1:(length(breaks)-1), function(i){
    win=(breaks[i]+1):breaks[i+1]
    id = unique(snphapid[win])
    if(length(id)!=1) stop("multiple id for a single window")
    whichsnp1=win[1]
    whichsnp2=win[length(win)]
    
    nfounder=length(strsplit(id,";")[[1]])
    
    data.frame(whichsnp1= whichsnp1, whichsnp2= whichsnp2, founder=id, nfounder = nfounder)
  }))
  
  
  haplotype$rilname=rilname
  
  haplotype$pos1=info$pos[haplotype$whichsnp1]
  haplotype$pos2=info$pos[haplotype$whichsnp2]
  haplotype$cM1=info$cM[haplotype$whichsnp1]
  haplotype$cM2=info$cM[haplotype$whichsnp2]
  
  
  haplotype
  
}


# STUFF I NEED TO WORK ON
if(F){
  
  #### FUSE close haplotypes Separated by only one SNP  
  
  
  
  #HAPLO = HAPLO[order(HAPLO$pos1,HAPLO$pos2),]
  
  fused = do.call(rbind, lapply(1:nrow(HAPLO), function(i){
    #print(i)
    whichfounder = which(HAPLO[i,3:ncol(HAPLO)]>0)
    founder = matrix(unlist(round(HAPLO[,whichfounder+2])), nrow = nrow(HAPLO))
    founder[founder>0.0] = 1
    nextcommon =  sort(unlist( apply(founder, 2, function(x){which(x>0)}) ))
    nextcommon = min(nextcommon[nextcommon >i])
    if(nextcommon==Inf){return(NULL)}
    
    j = nextcommon 
    
    overlapwithnext = HAPLO[j,1] <= HAPLO[i,2]
    
    whichfounder = HAPLO[i,3:ncol(HAPLO)]
    whichfounder2 = HAPLO[j,3:ncol(HAPLO)]
    whichfounder[whichfounder>0] = round(1, digits=0)
    whichfounder2[whichfounder2>0] = round(1, digits=0)
    
    commonfounder = which((whichfounder + whichfounder2) > 1)
    
    
    if(!overlapwithnext){
      #print(i)
      
      s1 = HAPLO[i,1]
      e1 = HAPLO[i,2]
      s2 =HAPLO[j,1]
      e2 =HAPLO[j,2]
      #print(s1)
      
      gap = e1:s2
      
      divergent = abs(founders[gap, commonfounder] - ril[gap])
      divergent = matrix(unlist(divergent), nrow = length(gap))
      divergent  = unique(apply(divergent, 2, function(x){which(x==1)}))
      
      gaprelsize = diff(snps$pos[gap[c(1,length(gap))]])/diff(snps$pos[c(b1,b2)])
      
      
      if(gaprelsize < 0.03 & length(divergent)<=1){
        win = c(s1:e1, s2:e2)
        
        new = runlsei(gx=ril[win],predictors=founders[win,])
        new = matrix(c(s1,e2,new), nrow=1)
        colnames(new)=c("pos1", "pos2", 1:ncol(founders))
        return(new)
      }else{return(NULL)}
      
    }else{return(NULL)}
    
  }))
  
  #if(is.null(fused)){fused = HAPLO[which(F),]}
  
  ### Fuse the fused haplotype that overlap
  i = 1
  while(i < nrow(fused)){
    j=i+1
    #pos1=HAPLO[i,1];pos2=HAPLO[i,2]
    overlapwithnext = fused[j,1] <= fused[i,2]
    
    x1 = fused[i,3:ncol(fused)]
    x2 = fused[j,3:ncol(fused)]
    x1[x1>0] = round(1, digits=0)
    x2[x2>0] = round(1, digits=0)
    
    commonfounder = x1 == x2
    
    if(sum(!commonfounder)==0 & overlapwithnext){
      b1 = fused[i,1]
      b2 =fused[j,2]
      x = fused[i,3:ncol(fused)]
      x = c(b1,b2,x)
      fused[i,] = x
      fused = fused[-j,]
      fused = matrix(fused, ncol=ncol(HAPLO))
    }else{i = i+1}
  }
  
  # Bind with other haplotypes & supress overlaping
  colnames(fused) = colnames(HAPLO)
  HAPLO = rbind(HAPLO, fused)
  HAPLO=supressEncased(HAPLO)
  
  
  
  #########################################
  ##### When two haplotype at the same place, chose 1
  
  i = 2
  while(i<nrow(HAPLO)){
    
    b1  = HAPLO[i,1]
    b2 = HAPLO[i,2]
    
    overlap1 = min(which(HAPLO[,1] < b1 & HAPLO[,2] > b1))
    overlap2 = max(which(HAPLO[,1] < b2 & HAPLO[,2] > b2))
    
    sameplace = which(HAPLO[,1] > HAPLO[overlap1,1] & HAPLO[,2] < HAPLO[overlap2,2])
    sameplace = sameplace[sameplace>i]
    
    if(length(sameplace)>0){
      keep = which.max(HAPLO[c(i,sameplace), 2] -HAPLO[c(i,sameplace), 1])
      out = c(i,sameplace)[which(!(1:(length(sameplace)+1) %in% keep))]
      HAPLO = HAPLO[-out,]
    }
    #print(i)
    #print(nrow(HAPLO))
    i=i+1
    
    
  }
  
  
  
  HAPLO = HAPLO[order(HAPLO$pos1,HAPLO$pos2),]
  HAPLO[,1] = snps$pos[HAPLO[,1]]
  HAPLO[,2] = snps$pos[HAPLO[,2]]
  
  for(i in 1:(nrow(HAPLO)-1)){
    j = i + 1
    sizei = HAPLO[i,2]- HAPLO[i,1]
    sizej = HAPLO[j,2] - HAPLO[j,1]
    
    wj = sizei/(sizei+sizej)
    wi = sizej/(sizei+sizej)
    
    newlimit = round(weighted.mean(c(HAPLO[i,2],HAPLO[j,1]), c(wi,wj)))
    
    HAPLO[i,2] = newlimit
    HAPLO[j,1] = newlimit -1
  }
}





