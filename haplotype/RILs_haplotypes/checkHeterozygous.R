library(readr)

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




getdis = function(win, line, founders){
  
  rx = line[win]
  fgt = founders[win,]
  
  # Get all the possible unique pairwise combinations of the targetted haplotypes
  # i.e. two haplotype within a founder
  comb = expand.grid(1:ncol(fgt), 1:ncol(fgt))
  comb = comb[comb$Var2 >= comb$Var1,]
  colnames(comb)=c("g1","g2")
  
  # For each pair, caluclate the genotype distance for each n snp with the poolseq infered ones
  # snpdist is a matrix where column are snp, rows corresponds to the possible combinations in combs
  # and value difference between diploid genotype
  minsnpdist=min(unlist(lapply(1:nrow(comb), function(i){
    g1 = comb[i,1] 
    g2 = comb[i,2] 
    cgt = (fgt[,g1]+fgt[,g2])/2 #
    
    #dist = sum(abs(cgt-fx), na.rm=T)
    sum(abs(cgt-rx),na.rm=T)
    
  })))
  
  edist = sum(rx[rx==0.5], na.rm=T)
  
  data.frame(whichsnp1 = min(win),whichsnp2 = max(win), genotypedistance=minsnpdist,
             nsnp=sum(!is.na(rx)),
             expectedDist = edist,
             phet = sum(rx==0.5,na.rm=T)/sum(!is.na(rx)))
  
}



getdis = function(win, line, founders){
  
  rx = line[win]
  fgt = founders[win,]
  
  # Get all the possible unique pairwise combinations of the targetted haplotypes
  # i.e. two haplotype within a founder
  comb = expand.grid(1:ncol(fgt), 1:ncol(fgt))
  comb = comb[comb$Var2 >= comb$Var1,]
  colnames(comb)=c("g1","g2")
  
  # For each pair, caluclate the genotype distance for each n snp with the poolseq infered ones
  # snpdist is a matrix where column are snp, rows corresponds to the possible combinations in combs
  # and value difference between diploid genotype
  minsnpdist=min(unlist(lapply(1:nrow(comb), function(i){
    g1 = comb[i,1] 
    g2 = comb[i,2] 
    cgt = (fgt[,g1]+fgt[,g2])/2 #
    
    #dist = sum(abs(cgt-fx), na.rm=T)
    sum(abs(cgt-rx),na.rm=T)
    
  })))
  
  data.frame(whichsnp1 = min(win),whichsnp2 = max(win), genotypedistance=minsnpdist)
  
}



getbadsnps = function(win, line, founders){
  
  
  rx = line[win]
  fgt = founders[win,]
  
  # Get all the possible unique pairwise combinations of the targetted haplotypes
  # i.e. two haplotype within a founder
  comb = expand.grid(1:ncol(fgt), 1:ncol(fgt))
  comb = comb[comb$Var2 >= comb$Var1,]
  colnames(comb)=c("g1","g2")
  
  
  snpdist=unlist(lapply(1:nrow(comb), function(i){
    g1 = comb[i,1] 
    g2 = comb[i,2] 
    cgt = (fgt[,g1]+fgt[,g2])/2 #
    
    #dist = sum(abs(cgt-fx), na.rm=T)
    sum(abs(cgt-rx),na.rm=T)
    
  }))
  
  min(snpdist)
  
  thiscomb = comb[which.min(snpdist),]
  g1=thiscomb[,1]
  g2=thiscomb[,2]
  cgt = (fgt[,g1]+fgt[,g2])/2
  badsnps = which(abs(cgt-rx)>0)+min(win)-1
  
  badsnps
  
}

setwd("/Users/tomparee/Documents/Documents - MacBook Pro de tom/")
load(file="./rockmanlab/becei/haplotype/founderHaplotypes/beceiFounderPhasedHaplotypes_chrI_complete.Rdata")
rils_unfiltered <- as.data.frame(read_csv("./rockmanlab/becei/genotype/Unfiltered/chrI_genotypeTable_beceiRILs_unfiltered.csv"))
rils_unfiltered = rils_unfiltered[match(phasedfounderhaplotypes$POS,rils_unfiltered$pos),]


#phet = apply(rils_unfiltered, 2, function(xx){ sum(xx==0.5)})/nrow(rils_unfiltered)
# sort(phet)


rilnames = colnames(rils_unfiltered)
rilnames= rilnames[grepl("A_|B_",rilnames)]
#rilnames=c("A_QG3393","A_QG3340", "A_QG3229","A_QG3339", "B_QG3375", "A_QG3130", "B_QG3162","B_QG3460", "B_QG3204", "A_QG3265")
hetgdis = do.call(rbind, lapply(rilnames, function(thisril){
  
  #thisril = "A_QG3340"
  print(thisril)
  line = rils_unfiltered[,thisril]
  linehom = line
  linehom[linehom==0.5] = NA
  windows = get.win(size=length(line), winsize=1000, minsize=1000)
  
  if(grepl("A_",thisril)) foundernames = c("FA.g1","FA.g2", "FM.g1", "FM.g2")
  if(grepl("B_",thisril)) foundernames = c("FB.g1","FB.g2", "FM.g1", "FM.g2")
  
  # print(foundernames)
  
  OUTPUT = do.call(rbind, lapply(windows, function(win){
    #win = windows[[12]]
    print(paste0("Progress: ",round((100*win[1])/length(line), digits=2), "%"))
    
    GD = getdis(win=win, line=line, founders=phasedfounderhaplotypes[,foundernames])
    
    if(GD$genotypedistance  > 0){
      
      ix =1:length(win)
      v =unlist(lapply(ix, function(splithere){
        
        w1 = win[1:splithere]
        w2 = win[splithere:length(win)]
        
        gd1 = getdis(win=w1, line=line, founders=phasedfounderhaplotypes[,foundernames])$genotypedistance
        gd2 = getdis(win=w2, line=line, founders=phasedfounderhaplotypes[,foundernames])$genotypedistance
        
        gd1+gd2
      }))
      
      if(min(v)< GD$genotypedistance){
        
        breakpoint = which.min(v)
        win1 = win[1:breakpoint]
        win2 = win[breakpoint:length(win)]
        
        GD  = rbind(getdis(win=win1, line=line, founders=phasedfounderhaplotypes[,foundernames]),
                    getdis(win=win2, line=line, founders=phasedfounderhaplotypes[,foundernames]))
      }
      
      
    }
    
    
    GD=cbind(GD,do.call(rbind,lapply(1:nrow(GD), function(i){
      w = GD$whichsnp1[i]:GD$whichsnp2[i]
      genotypedistanceHom = getdis(win=w, line=linehom, founders=phasedfounderhaplotypes[,foundernames])$genotypedistance
      nhet = sum(line[w]==0.5)
      nsnp = sum(!is.na(line[w]))
      
      #badsnpshom = getbadsnps(win=w, line=linehom, founders=phasedfounderhaplotypes[,foundernames])
      #badsnpshet = getbadsnps(win=w, line=line, founders=phasedfounderhaplotypes[,foundernames])
      
      #paste(badsnpshom,collapse=';')
      #paste(badsnpshet,collapse=';')
      
      getdis(win=w, line=line, founders=phasedfounderhaplotypes[,foundernames])
      
      data.frame(genotypedistanceHom,nhet,nsnp)
      
    })))
    
    
    
    return(GD)
    
    
  }))
  
  OUTPUT$ril = thisril
  
  OUTPUT
  
  
}))



test = subset(hetgdis, ril == rilnames[[10]])


ggplot(test)+
  geom_step(aes(x=whichsnp1, y=nhet*0.5), color='orange',size=1.7, alpha=0.5)+
  geom_segment(aes(x=whichsnp1,xend = whichsnp2, y=nhet*0.5, yend=nhet*0.5), color='orange', size=1/7, alpha=0.5)+
  geom_step(aes(x=whichsnp1, y=genotypedistance))+
  geom_segment(aes(x=whichsnp1,xend = whichsnp2, y=genotypedistance, yend=genotypedistance ))+
  geom_step(aes(x=whichsnp1, y=genotypedistanceHom), color='darkgreen')+
  geom_segment(aes(x=whichsnp1,xend = whichsnp2, y=genotypedistanceHom, yend=genotypedistanceHom), color='darkgreen')+
  xlab("Position (ordered snps)")+
  ylab("Genotype distance")+
  ggtitle(paste0(test$ril[1],"; het:", round(phet[test$ril[1]]*100, digits = 2), "%"))#+
  #coord_cartesian(ylim=c(0,10))

gridExtra::grid.arrange(p1,p2,p3,nrow=3)


