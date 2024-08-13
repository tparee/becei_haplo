library(readr)

#########################################################################
################ FUNCTIONS ##############################################


load(file="~/Documents/Documents - MacBook Pro de tom/rockmanlab/cbecei_haplo/rilsAfounderHaplotypeBlocks_chrI.Rdata")
RILsFounderHaploBlocks = haplotype
  
countBreaks = function(windows, rils, foundergt,info){
  do.call(rbind, lapply(1:(length(windows)-1), function(i){
    
    win = windows[i]:windows[i+1]
    if(i %% 1000 ==0 ) print(i)
    fx = foundergt[win,]
    gx = rils[win,]
    
    nbreak = sum(apply(gx, 2, function(x){ ifelse(sum(apply(fx!=x,2,mean,na.rm=T) == 0) == 0, 1,0) }), na.rm=T)
    data.frame(nbreak = nbreak, gsize = info$cM[win[length(win)]] - info$cM[win[1]], start = win[1], end = win[length(win)])
  }))
}


whichRilsBroken = function(rilsgt, foundergt){
  which(apply(rilsgt, 2, function(x){ ifelse(sum(apply(foundergt!=x,2,mean,na.rm=T) == 0) == 0, T,F) }))
  }




#rils = cbind(rilsA,rilsB)
Search.LowLDSNP = function(rils, winsizeLD = 500, LDth=0.9, min.nHighlink=3){
  
  steps = seq(1, nrow(rils), winsizeLD/2)
  steps = steps[-1]
  
  whichlowLDSNPs = unlist(lapply(steps, function(step){
    
    print(step)
    win = (step-winsizeLD):(step+winsizeLD)
    win = win[win >= 1 & win <= nrow(rils)]
    rx = cor(t(rils[win,]), use = "pairwise.complete.obs")
    
    #For each SNP count the number of other SNP it is in LD with
    nHighlink = apply(rx^2, 2, function(x){sum(x>LDth, na.rm=T)}) -1 #-1 because to account for the diagonal
    
    win[which(nHighlink<min.nHighlink)]
  }))
  
  whichlowLDSNPs=sort(unique(whichlowLDSNPs))
  whichlowLDSNPs
}







#InferFoundersHaploBlocksWIN(win,founderA, founderB, founderM, rilsA, rilsB, info)

InferFoundersHaploBlocksWIN = function(win,founderA, founderB, founderM, rilsA, rilsB, info){
  #fgt = founder genotypes infered by poolseq within the target window (1=hom alt; 0.5=het; 0 = hom ref)
  fgt = list(A=(founderA[win,1]+founderA[win,2])/2, 
             B=(founderB[win,1]+founderB[win,2])/2,
             M= (founderM[win,1]+founderM[win,2])/2)
  
  # List crontaining rils genotype within the window for cross A and cross B
  rilgeno = list(A=rilsA[win,], B=rilsB[win,])
  
  # Infer founder haplotype in each cross (=major haplotypes in rils)
  rilshap = lapply(rilgeno, function(gx) {
    
    # Remove RILs with a high proportion of missing values
    gx = gx[, apply(gx, 2, function(x) { sum(is.na(x)) }) / nrow(gx) < 0.3]
    
    # Identify major haplotypes in rils
    grouprils = cutree(hclust(dist(t(gx))), h = 0) # Cluster RILs based on their haplotypes
    npergroup = table(grouprils) # Count the number of RILs in each cluster
    npergroup = npergroup[npergroup > 10] # Keep only clusters with more than 10 RILs (I don't have a good rational for this number, could be more)
    npergroup = sort(npergroup, decreasing = TRUE)  # Sort the clusters by size in decreasing order
    if (length(npergroup) > 4) npergroup = npergroup[1:4] # Limit to the top 4 major haplotypes, if more than 4 are present (we expect four founder haplo at max)
    
    # Gat each major haplotype 
    # (Use info from all RILs from a cluster rather than one, though they are identical, because they are some NA)
    rilshapx = do.call(cbind, lapply(as.numeric(names(npergroup)), function(x) {
      apply(gx[, grouprils == x], 1, mean, na.rm = TRUE)
    }))
    
    rilshapx
  })
  
  
  
  
  # First, let's find the common RILs haplotypes between cross A and cross B (should correspond to the founder M)
  # Here we generate a df with hA the ith column number in cross A, hB the ith column number in cross B, and percentMatch the % of match between the two (0=>1)
  sharedRilsHaplo = setNames(as.data.frame(t(apply(expand.grid(1:ncol(rilshap[["A"]]), 1:ncol(rilshap[["B"]])), 1, function(pair){
    hA = rilshap[["A"]][,pair[1]]
    hB = rilshap[["B"]][,pair[2]]
    percentMatch = sum(hA==hB, na.rm=T)/sum(!is.na(hA) & !is.na(hB))
    c(pair,percentMatch)
  }))), c("hA", "hB", "percentMatch"))
  # Only keep pairs of haplotypes with percentMatch == 1, i.e, haplotype that are the same between cross A and B
  sharedRilsHaplo = sharedRilsHaplo[sharedRilsHaplo$percentMatch==1,]
  
  
  # Classify the RILs haplotype in shared and unshared: 
  rilshap = list(onlyA=matrix(rilshap[["A"]][,-sharedRilsHaplo$hA], nrow=nrow(rilshap[["A"]])), #Rils haplotypes only found in cross A
                 onlyB=matrix(rilshap[["B"]][,-sharedRilsHaplo$hB], nrow=nrow(rilshap[["B"]])), #Rils haplotypes only found in cross B
                 shared = matrix(rilshap[["B"]][,sharedRilsHaplo$hB], nrow=nrow(rilshap[["B"]]))) # Shared Rils haplotypes
  
  # Transform in a matrix with column names indicating shared, onlyA or onlyB
  rilshap = do.call(cbind, lapply(1:length(rilshap), function(i){
    x=rilshap[[i]]
    if(ncol(x)>0) colnames(x) = paste0(names(rilshap)[i], 1:ncol(x)) 
    x
  }))
  
  
  
  # For each founder, calulcate the genotype distance (dist) 
  # Between the genotype infered from poolseq (stored within fgt) and the possible pairs of haplotypes
  # We focus on the haplotype that are relevant for each founder (i.e, founder A, we only look at haplo identified in cross A or shared)
  # A result for each founder in returned into a list with genotype ranked  by genetic distances
  gdistance = lapply(c("A", "B", "M"), function(thisfounder){
    
    fx = fgt[[thisfounder]]
    
    # In rilsHhap, restrict the haplotypes that corresponds to this founder
    # If founder M = shared haplo; if A/B, only cross A/B or shared
    targetHaploRils = c("onlyA", "onlyB", "shared")[which(c("A", "B", "M")==thisfounder)]
    targetHaploRils = which(grepl(paste0(targetHaploRils,"|shared"), colnames(rilshap)))
    
    # Get all the possible unique pairwise combinations of the targetted haplotypes
    # i.e. two haplotype within a founder
    comb = expand.grid(targetHaploRils, targetHaploRils)
    comb = comb[comb$Var2 >= comb$Var1,]
    colnames(comb)=c("g1","g2")
    
    # For each pair, caluclate the genotype distance for each n snp
    # i.e., difference between diploid genotype 
    snpdist=do.call(rbind,(lapply(1:nrow(comb), function(i){
      g1 = comb[i,1] 
      g2 = comb[i,2] 
      cgt = (rilshap[,g1]+rilshap[,g2])/2
      
      #dist = sum(abs(cgt-fx), na.rm=T)
      matrix(abs(cgt-fx), nrow=1)
    })))
    
    # get the genotype distance summed across all sites within the window for each combinations
    windist = apply(snpdist, 1, sum, na.rm=T)
    
    comb = cbind(comb,windist = windist , snpdist)
    colnames(comb)[4:ncol(comb)] = paste0("snpdist", colnames(comb)[4:ncol(comb)])
    #Order by increasing distance
    # then by ifelse(comb[,"g2"]-comb[,"g1"] ==0,0,1) => This just ensure that if different combinations have the same distance,
    # The combinations with different rils haplotype will be classed higher => maximize the number of rils haplotype represented
    comb = comb[order(comb[,"windist"], ifelse(comb[,"g2"]-comb[,"g1"] ==0,0,1), decreasing = c(F,T)),]
    
    # matrix with colum = g1,g2: the rils haplotype combination (value = ith column in rilsname)
    # windist the genotype distance summed across the windo
    #snpsdist1,2,.. = snp distance for each individual snp
    comb
    
  })
  
  # Now we want to choose the haplotype attribution that minimize the total genetic distance 
  # i.e. taking the first combination (already ranked by distance) for each founder
  
  # Now we want to choose the haplotype attribution that minimize the total genetic distance 
  # i.e. taking the first combination (already ranked by distance) for each founder
  
  #gdistance2 = lapply(gdistance, function(x){x[,1:3]})
  fhap = unlist(lapply(gdistance, function(x){x[1,c(1,2)]}))
  totdist = sum(unlist(lapply(gdistance, function(x){x[1,3]})))
  
  do.call(rbind,lapply(gdistance, function(x){x[1,4:ncol(x)]}))

  names(fhap) = c('FA.g1', 'FA.g2', 'FB.g1', 'FB.g2', 'FM.g1', 'FM.g2')
  
  # Verify that all the major haplotypes identified within RILs are attributed in 
  representedRilsHap = sort(unique(fhap))
  allRepresented = sum(1:ncol(rilshap) %in% representedRilsHap)==ncol(rilshap)
  
  
  if(!allRepresented){
    #print(step)
    
    #The non represented haplotype(s)
    missingHaplo = which(is.na(match(1:ncol(rilshap), representedRilsHap)))
    
    # Generate a table with the number of different snp between the missing haplotypes and other (non-missing) haplotypes
    hapDiff = setNames(as.data.frame(t(apply(expand.grid(missingHaplo, (1:ncol(rilshap))[-missingHaplo]), 1, function(pair){
      h1 = rilshap[,pair[1]]
      h2 = rilshap[,pair[2]]
      ndiff = sum(h1!=h2, na.rm=T)
      c(pair,ndiff)
    }))), c("missHap", "otherHap", "ndifferent"))
    
    # If only one different SNP between the missing hap and another hap,
    # consider that the two haplotype are the same and that the divergent snp is problematic and exclude it
    ## Note: Here we try to infer a robust base for the founder haplotypes with only SNP that we are very confident with
    ## there will be time later to add low quality SNP and other suspicious SNP, so SNP exluded here are not lost forever.
    
    hapDiff = hapDiff[hapDiff$ndifferent==1,]
    
    if(nrow(hapDiff)>0){
      badsnp = unlist(lapply(1:nrow(hapDiff), function(i){
        x=hapDiff[i,]
        h1 = rilshap[,unlist(x[1])]
        h2 = rilshap[,unlist(x[2])]
        which(h1!=h2)
      }))
      
      if(nrow(hapDiff) != length(badsnp)){print("DEBUG ME")}
      
      allRepresented = T
      win = win[-badsnp]
      rilshap = rilshap[-badsnp,]
      rilgeno = lapply(rilgeno, function(x){x[-badsnp,]})
      totdist=sum(unlist(lapply(gdistance, function(x){x[1,-(3+badsnp)]})))
    }
    
  }
  
  whichrilshap = fhap
  het = ifelse(fhap[c(1,3,5)]-fhap[c(2,4,6)] == 0,F,T)
  fhap = rilshap[,fhap]
  colnames(fhap) = c('FA.g1', 'FA.g2', 'FB.g1', 'FB.g2', 'FM.g1', 'FM.g2')
  rilgeno=do.call(cbind, lapply(rilgeno, function(x){x}))
  nbreaks=countBreaks(windows=c(1,length(win)), rilgeno, fhap, info[win,])$nbreak
  #nbreaks2=sum(!(grouprils %in% (as.numeric(names(npergroup)))))
  #if(nbreaks1 > nbreaks2){print("PROBLEM");print(step)}
  
  return(list(win=win,foundergt=fhap, heterozygous = het, rilshap=whichrilshap, nbreaks = nbreaks, totdist = totdist, allRepresented=allRepresented, gsize=diff(range(info$cM[win]))))
}



i=647
win = haplotype[[i]]$win
wb = whichRilsBroken(rils[win,], fhap)

gx=rils[win,wb]
grouprils = cutree(hclust(dist(t(gx))), h = 0) # Cluster RILs based on their haplotypes
table(grouprils)

View(cor(t(gx), use="pairwise.complete.obs"))


extendedwin = unlist(lapply(haplotype[c((i-2):(i-1),(i+1):(i+2))], function(x){x$win}))


gx=rils[extendedwin,wb]
grouprils = cutree(hclust(dist(t(gx))), h = 0) # Cluster RILs based on their haplotypes
table(grouprils)






countBreaks(windows=c(1:length(win)), rilgeno, fhap, info[win,])$nbreak



InferFoundersHaploBlocks = function(founderA, founderB, founderM, rilsA, rilsB, info,nsnpWin = 50, maxSizeCM=3){
  
  
  # Define a sequence of steps for iterating over the RILs dataset based on a window size
  steps = seq(1, nrow(info), nsnpWin)
  steps[length(steps)] = nrow(info) + 1 
  
  # Ensure that genetic distance between adjacent markers is not too large
  # If recombination too high, it is more difficult to infer founder haplotypes because they are broken
  steps = c(unlist(lapply(1:(length(steps) - 1), function(step) {
    
    win = steps[step]:(steps[step + 1] - 1) # Define the window of markers for this step
    gsize = diff(range(info$cM[win])) # Calculate the genetic size of the window
    
    # Check if the genetic size is under the limit
    if (gsize < maxSizeCM) {
      # If within limits, return the first marker index of the window
      return(win[1])
    } else {
      # If not, return the indices of the first and last markers that maximize the genetic distance
      return(win[c(1, which.max(diff(info$cM[win])) + 1)])
    }
  })), steps[length(steps)])
  
  
  

  founderhaplotypes = lapply(1:(length(steps)-1), function(step){
    #print(step)
    if(step %in% 10 == 0) print(step)
    
    # The window under consideration (this step to the next one -1)
    win =steps[step]:(steps[step+1]-1)
    
    # Infer founder haplotype blocks in the window + some stat about inference
    InferFoundersHaploBlocksWIN(win=win,
                                founderA=founderA, founderB=founderB, founderM=founderM,
                                rilsA=rilsA, rilsB=rilsB, 
                                info)
  })
  
  
  
  return(founderhaplotypes)
}

which(unlist(lapply(haplotype, function(x){x$win[1]}))>=30701)
which(unlist(lapply(founderhaplotypessave, function(x){x$win[1]}))>=30701)

#founderhaplotypes=founderhaplotypessave
#founderhaplotypes = haplotype
phaseHaplotypes = function(founderhaplotypes, info, founderA, founderB, founderM){
  
  i = 1
  #while(i < 379){
  while(i < length(founderhaplotypes)-1 ){
    print(i)
    
    fi = founderhaplotypes[[i]]$foundergt
    fj = founderhaplotypes[[i+1]]$foundergt
    wini = founderhaplotypes[[i]]$win
    winj = founderhaplotypes[[i+1]]$win
    midwin = c(wini[floor(length(wini)/2):length(wini)], winj[1:floor(length(winj)/2)])
    winall = c(wini, winj)
    #winall = c(winall, max(winall)+1) # Just some add on to compensate a dumb thing that I am to lazy to correct nicely right now
    
    potentialSwap = setNames(expand.grid(c(1,-1), c(1,-1), c(1,-1)), c("swapFA", "swapFB", "swapFM"))
    
    nbreaks = unlist(lapply(1:nrow(potentialSwap), function(ii){
      
      if(potentialSwap[ii,1]==1){swapFA=c(1,2)}else{swapFA=c(2,1)}
      if(potentialSwap[ii,2]==1){swapFB=c(3,4)}else{swapFB=c(4,3)}
      if(potentialSwap[ii,3]==1){swapFM=c(5,6)}else{swapFM=c(6,5)}
      swap = c(swapFA,swapFB, swapFM)
      
      #fused = rbind(fi[wini %in% midwin,], fj[winj %in% midwin,swap])
      #countBreaks(windows=c(1,nrow(fused)), rils=rils[midwin,], foundergt=fused, info=info[midwin,])$nbreak
      
      fused = rbind(fi, fj[,swap])
      countBreaks(windows=c(1,nrow(fused)), rils=rils[winall,], foundergt=fused, info=info[winall,])$nbreak
    }))
    
    swap = as.data.frame(potentialSwap[which.min(nbreaks),])
    swap = c(apply(swap, 2, function(x){if(x==1){c(1,2)}else{c(2,1)}})) + c(0,0,2,2,4,4)
    
    haploall_phased = rbind(fi, fj[,swap])
    #nbreaksPhased = nbreaks[which.min(nbreaks)]
    nbreaksPhased = countBreaks(windows=c(1,nrow(haploall_phased)), rils=rils[winall,], foundergt=haploall_phased, info=info[winall,])$nbreak
    
    haploall_infered = InferFoundersHaploBlocksWIN(win=winall,
                                                   founderA=founderA,
                                                   founderB=founderB,
                                                   founderM=founderM,
                                                   rilsA=rilsA,
                                                   rilsB=rilsB,
                                                   info)
    
    #countBreaks(windows=c(1:nrow(haploall_phased)), rils=rils[winall,], foundergt=haploall_phased, info=info[winall,])$nbreak
    
    
      
    nbreaksInfered = haploall_infered$nbreaks
    
    haploall_infered = haploall_infered$foundergt
    
    choosePhased =  nbreaksInfered >= nbreaksPhased
    
    if(choosePhased){
      founderhaplotypes[[i+1]]$foundergt = fj[,swap]
      i=i+1
    }else{
      founderhaplotypes[[i+1]]$foundergt = haploall_infered
      founderhaplotypes[[i+1]]$win =  c(wini, winj)
      founderhaplotypes[[i+1]]$nbreaks =  nbreaksInfered
      founderhaplotypes = founderhaplotypes[-i]
      i=i-1
    }
    
    
  }
  
  #founderhaplotypes[[380]]
  whichsnp = unlist(lapply(founderhaplotypes, function(x){x$win}))
  phased = do.call(rbind, lapply(founderhaplotypes, function(x){x$foundergt}))
  
  founderPoolgt = cbind(founderA, founderB, founderM)[whichsnp,]
  founderPoolgt = (founderPoolgt[,c(1,3,5)] + founderPoolgt[,c(2,4,6)])/2
  snpdist = ((phased[,c(1,3,5)] + phased[,c(2,4,6)])/2)-founderPoolgt
  snpdist = apply(snpdist, 1, function(x){abs(sum(x))})
  phased = cbind(info[whichsnp,],Poolgenotype.dist=snpdist, phased)
  
  return(phased)
  
  
}









test = do.call(rbind, lapply(founderhaplotypes, function(x){data.frame(nbreaks=x$nbreaks,
                                                                totdist = x$totdist,
                                                                gsize = x$gsize,
                                                                winstart = x$win[1])}))


test2 = do.call(rbind, lapply(haplotype, function(x){data.frame(nbreaks=x$nbreaks,totdist = x$totdist,gsize = x$gsize, winstart = x$win[1])}))


ggplot(test, aes(gsize, nbreaks))+geom_point()

#########################################################################
#########################################################################
rils <- read_csv("~/Documents/Documents - MacBook Pro de tom/rockmanlab/cbecei_haplo/041724_becei_chr_I_Cross_AB_Founder_RIL_GT_table.csv")
#rilsA <- read_csv("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/cbecei_haplo/becei_chr_I_Cross_A_filtered_genotype_table.csv")
#rilsB <- read_csv("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/cbecei_haplo/becei_chr_I_Cross_B_filtered_genotype_table.csv")

founders = rils[,c('FM', 'FA', 'FB')]
info = rils[,1:4]

rils = rils[,grepl("_QG", colnames(rils))]
rilnames =  colnames(rils)



# some formating

#rils = rils[,7:ncol(rils)]
rils[rils=="1/1"] = "1"
rils[rils=="0/0"] = "0"
rils[rils=="0/1"]=NA
rils = matrix(as.numeric(unlist(c(rils))), ncol=ncol(rils))
colnames(rils) = rilnames


# remove fixed snps
fixed = apply(rils, 1, function(x){x = mean(x, na.rm=T); x %in% c(0,1)})
nafreq = apply(rils,1, function(x){sum(is.na(x))})/ncol(rils)

rils = rils[!fixed & nafreq<0.05,]
info = info[!fixed & nafreq<0.05,]
founders = founders[!fixed & nafreq<0.05,]

whichlowLDSNPs = Search.LowLDSNP(rils, winsizeLD = 500, LDth=0.9, min.nHighlink=3)
rils = rils[-whichlowLDSNPs,]
info = info[-whichlowLDSNPs,]
founders = founders[-whichlowLDSNPs,]


#badsnps = "I:3851354:T-G" #snp that is really suspicious (lot of NA value and cause many breakpoints which is consistent with genotyping error)
#rils[(info$ID %in% badsnps), ]=NA
#founders[(info$ID %in% badsnps), ] = NA

rilsA = rils[,grepl("A_QG", colnames(rils))]
rilsB = rils[,grepl("B_QG", colnames(rils))]


#founder 1 & 2 in format snp x two haplotype
founderM = founders$FM
founderA = founders$FA
founderB = founders$FB
founderM = do.call(rbind, lapply(strsplit(founderM, "/"), function(x){as.numeric(x)}))
founderA = do.call(rbind, lapply(strsplit(founderA, "/"), function(x){as.numeric(x)}))
founderB = do.call(rbind, lapply(strsplit(founderB, "/"), function(x){as.numeric(x)}))


# First, InferFoundersHaploBlocks look at the different possible haplotypes in each genomic window
# the most frequent haplotypes are assumed to be acestral
# each window need to be big enough to diferenciate founders  (default is 50SNP)
# And not too large in genetid distance => if too much recombination, the ancestral haplotype are too broken to be recognize
# If it is the case, the window is subdivided
# Then the major haplotypes are attributed to each founders in the way that minimize the change from the initially inferred genotypes (given in founder1 and founder2)
# The function return the four founder hap in  c("founder1.g1","founder1.g2","founder2.g1","founder2.g2") order for each window, in a list
founderhaplotypes = InferFoundersHaploBlocks(founderA, founderB, founderM, rilsA, rilsB, info,nsnpWin = 100, maxSizeCM=3)
#founderhaplotypessave = founderhaplotypes
#This function phase the haplotype block inferred above in a way that minimize the number of breakpoints
phasedfounderhaplotypes = phaseHaplotypes(founderhaplotypes, info, founderA, founderB, founderM)

phasedfounderhaplotypes = phasedfounderhaplotypes[-(30752:30762),]

#Problematic snp for chromosome I: 30752:30762

phasedfounderhaplotypes=cbind(info[phasedfounderhaplotypes$whichsnp,], phasedfounderhaplotypes[,-1])

save(phasedfounderhaplotypes, file="~/Documents/Documents - MacBook Pro de tom/rockmanlab/cbecei_haplo/240501_beceiFounderPhasedHaplotypes_chrI.Rdata")



#As a "control" we can count the number of breaks that happened in different intervals
#(= the haplotype that are not identical to one founder for a given window)
# window by snps:
#founderhaplotypes=phased

mm = match(founderhaplotypes$POS, info$POS)

windows = seq(1,nrow(founderhaplotypes), 100)
windows[length(windows)]=nrow(founderhaplotypes)
breakBins = countBreaks(windows, rils[mm,], founderhaplotypes[,5:10], info[mm,])
sum(breakBins$nbreak)
ggplot(breakBins, aes(gsize, nbreak))+geom_point()
# we are very close from the number of recombination event one may expect from 5 generation of oucrossing with 1 CO per meiosis
# expected: 148 lines * 5 generation * 0.5 CO per chrom = 370
# observed: 374










rils2 = cbind(rilsA, rilsB)[mm,]





win = 19000:22000

test = InferFoundersHaploBlocksWIN(adjwin,founderA[mm,], founderB[mm,], founderM[mm,], rilsA[mm,], rilsB[mm,], info[mm,])
test = test$foundergt

diff(range(founderhaplotypes$cM[win]))

adjwin = (min(win)-1000):(max(win)+1000)
adjwin = adjwin[!(adjwin %in% win)]

fhap = founderhaplotypes[adjwin,grepl("F",colnames(founderhaplotypes))]
rhap = rils2[adjwin,]

groups = cutree(hclust(dist(t( cbind(fhap, rhap) ))), h = 0) # Cluster RILs based on their haplotypes

fgroups = groups[1:ncol(fhap)]

table(groups)

View(cbind(rhap[,"B_QG3492"],fhap))

lapply(fgroups, function(thisfoundergroup){
  print()
  unbroken = which(groups %in% thisfoundergroup)
  unbroken = unbroken - ncol(fhap)
  unbroken=unbroken[unbroken>=1] 
  
  midhap = rils2[win,unbroken]
  midhapgrp = cutree(hclust(dist(t( midhap ))), h = 0)
  
  paste(c(sort(table(midhapgrp),decreasing=T)), collapse=";")
  
  
})


npergroup = table(grouprils) # Count the number of RILs in each cluster
npergroup = npergroup[npergroup > 10] # Keep only clusters with more than 10 RILs (I don't have a good rational for this number, could be more)
npergroup = sort(npergroup, decreasing = TRUE)  # Sort the clusters by size in decreasing order
if (length(npergroup) > 4) npergroup = npergroup[1:4]







win = 33601:33701

corfound = cor(t(founderhaplotypes[win,5:10]), use = "pairwise.complete.obs")
corrils = cor(t(rils2[win,]), use = "pairwise.complete.obs")

cor_ratio = corrils/corfound


apply(cor_ratio, 2, mean)


save(founderhaplotypes, file="~/Documents/Documents - MacBook Pro de tom/rockmanlab/cbecei_haplo/founderPhasedHaplotypes_chrI.Rdata")





#founderhaplotypes=founderhaplotypes2
#founderhaplotypes = founderhaplotypes[c(1:(length(founderhaplotypes)-1) %% 2 == 1, T)]

phaseHaplotypes = function(founderhaplotypes, info, founderA, founderB, founderM){
  
  i = 1
  #while(i < 379){
  while(i < length(founderhaplotypes)-1 ){
    print(i)
    
    fi = founderhaplotypes[[i]]$foundergt
    fj = founderhaplotypes[[i+1]]$foundergt
      
    wini = founderhaplotypes[[i]]$win
    winj = founderhaplotypes[[i+1]]$win
    #midwin = c(wini[floor(length(wini)/2):length(wini)], winj[1:floor(length(winj)/2)])
    winall = c(wini, winj)
    #winall = c(winall, max(winall)+1) # Just some add on to compensate a dumb thing that I am to lazy to correct nicely right now
    
    potentialSwap = setNames(expand.grid(c(1,-1), c(1,-1), c(1,-1)), c("swapFA", "swapFB", "swapFM"))
    
    nbreaks = unlist(lapply(1:nrow(potentialSwap), function(ii){
      
      if(potentialSwap[ii,1]==1){swapFA=c(1,2)}else{swapFA=c(2,1)}
      if(potentialSwap[ii,2]==1){swapFB=c(3,4)}else{swapFB=c(4,3)}
      if(potentialSwap[ii,3]==1){swapFM=c(5,6)}else{swapFM=c(6,5)}
      swap = c(swapFA,swapFB, swapFM)
      
      #fused = rbind(fi[wini %in% midwin,], fj[winj %in% midwin,swap])
      #countBreaks(windows=c(1,nrow(fused)), rils=rils[midwin,], foundergt=fused, info=info[midwin,])$nbreak
      
      fused = rbind(fi, fj[,swap])
      countBreaks(windows=c(1,nrow(fused)), rils=rils[winall,], foundergt=fused, info=info[winall,])$nbreak
    }))
    
    if(which.min(nbreaks)!=1) print("HERE")
    swap = as.data.frame(potentialSwap[which.min(nbreaks),])
    swap = c(apply(swap, 2, function(x){if(x==1){c(1,2)}else{c(2,1)}})) + c(0,0,2,2,4,4)
    
    haploall_phased = rbind(fi, fj[,swap])
    #nbreaksPhased = nbreaks[which.min(nbreaks)]
    nbreaksPhased = countBreaks(windows=c(1,nrow(haploall_phased)), rils=rils[winall,], foundergt=haploall_phased, info=info[winall,])$nbreak
    
    haploall_infered = InferFoundersHaploBlocksWIN(win=winall,
                                                   founderA=founderA,
                                                   founderB=founderB,
                                                   founderM=founderM,
                                                   rilsA=rilsA,
                                                   rilsB=rilsB,
                                                   info)
    
    #countBreaks(windows=c(1:nrow(haploall_phased)), rils=rils[winall,], foundergt=haploall_phased, info=info[winall,])$nbreak
    
    
    
    nbreaksInfered = haploall_infered$nbreaks
    
    haploall_infered = haploall_infered$foundergt
    
    choosePhased =  nbreaksInfered >= nbreaksPhased
    
    
    
    if(choosePhased){
      founderhaplotypes[[i+1]]$foundergt = fj[,swap]
      i=i+1
    }else{
      founderhaplotypes[[i+1]]$foundergt = haploall_infered
      founderhaplotypes[[i+1]]$win =  c(wini, winj)
      founderhaplotypes[[i+1]]$nbreaks =  nbreaksInfered
      founderhaplotypes = founderhaplotypes[-i]
      i=i-1
    }
    
    
  }
  
  #founderhaplotypes[[380]]
  whichsnp = unlist(lapply(founderhaplotypes, function(x){x$win}))
  phased = do.call(rbind, lapply(founderhaplotypes, function(x){x$foundergt}))
  
  founderPoolgt = cbind(founderA, founderB, founderM)[whichsnp,]
  founderPoolgt = (founderPoolgt[,c(1,3,5)] + founderPoolgt[,c(2,4,6)])/2
  snpdist = ((phased[,c(1,3,5)] + phased[,c(2,4,6)])/2)-founderPoolgt
  snpdist = apply(snpdist, 1, function(x){abs(sum(x))})
  phased = cbind(info[whichsnp,],Poolgenotype.dist=snpdist, phased)
  
  return(phased)
  
  
}


# # Extract the "statistics" for each haplo and store it in a data.frame
# ## (nbreaks,number of broken haplotype within RILs,
# ## totdist = sum of genetic distance between founders genotype infered by poolseq and their inferred haplo,
# ## allRepresented = T/F if all the major haplotypes found in the RILs could be attributed or not to the founders)
# hapstat = do.call(rbind, lapply(haplotype, function(x){data.frame(nbreaks=x$nbreaks,
#                                                                   totdist = x$totdist,
#                                                                   allRepresented=x$allRepresented,
#                                                                   gsize=x$gsize)}))
# 
# ggplot(hapstat, aes(gsize, nbreaks))+geom_point()
# # Suspicious haplotype that we going to examinate
# # We're going to be really strict and examine all haplotype that have a break in rils
# # or difference between founder genotype inferred from poolseq or rils haplotypes
# # or some where some major haplotypes are not attributed to one founder
# suspicious = which(hapstat$nbreaks >= 1 | hapstat$totdist > 0 | !hapstat$allRepresented)
# 
# 
# 
# suspicious = do.call(rbind, lapply(suspicious, function(i){
#   #i=89
#   
#   # Does extending the window "improve" the haplotype inference?
#   # Value can be changed, decision at the end
#   extend = F 
#   
#   extendedwin = unlist(lapply(haplotype[c((i-1):(i+1))], function(x){x$win})) # => window containing the adjacent haplotypes, including the suspicious haplotype
#   
#   #Haplo infered from extended windo
#   extHap = InferFoundersHaploBlocksWIN(win=extendedwin,
#                                        founderA=founderA, founderB=founderB, founderM=founderM,
#                                        rilsA=rilsA, rilsB=rilsB,
#                                        info)
#   
#   susp.nbreaks = hapstat$nbreaks[i]
#   ext.nbreaks = extHap$nbreaks
#   
#   if( ext.nbreaks < susp.nbreaks) extend = T
#   
#   data.frame(whichHaplo = i, extend = extend)
#   
#   
#   #adjacentwin =  unlist(lapply(haplotype[c((i-1),(i+1))], function(x){x$win})) # => window containing the adjacent haplotypes, but not the suspicious haplotype
#   
#   #suspiciousHap = InferFoundersHaploBlocksWIN(win=haplotype[[i]]$win,
#   #                                     founderA=founderA, founderB=founderB, founderM=founderM,
#   #                                     rilsA=rilsA, rilsB=rilsB,
#   #                                     info)
#   
#   #Haplo infered from adjacent window (exclude the susupicious haplotype)
#   #adjHap = InferFoundersHaploBlocksWIN(win=adjacentwin,
#   #                                     founderA=founderA, founderB=founderB, founderM=founderM,
#   #                                     rilsA=rilsA, rilsB=rilsB,
#   #                                     info)
#   
#   # adj.broken = whichRilsBroken(cbind(rilsA, rilsB)[adjHap$win,], adjHap$foundergt)
#   # ext.broken = whichRilsBroken(cbind(rilsA, rilsB)[extHap$win,], extHap$foundergt)
#   # 
#   # diffBroken = ext.broken[!(ext.broken %in% adj.broken)]
#   # 
#   # 
#   # diffBroken=11
#   # leftwin = extHap$win < min(extHap$win[!(extHap$win %in% adjHap$win)])
#   # rightwin = extHap$win > max(extHap$win[!(extHap$win %in% adjHap$win)])
#   # 
#   # left.gx = cbind(extHap$foundergt[leftwin,], cbind(rilsA, rilsB)[extHap$win[leftwin],diffBroken])
#   # right.gx = cbind(extHap$foundergt[rightwin,], cbind(rilsA, rilsB)[extHap$win[rightwin],diffBroken])
#   # adj.gx = cbind(extHap$foundergt[rightwin|leftwin,], cbind(rilsA, rilsB)[extHap$win[rightwin|leftwin],diffBroken])
#   # 
#   # 
#   # #ext.gx = cbind(extHap$foundergt, cbind(rilsA, rilsB)[extHap$win,diffBroken])
#   # 
#   # cutree(hclust(dist(t(adj.gx))), h = 0)
#   # cutree(hclust(dist(t(left.gx))), h = 0) # Cluster RILs based on their haplotypes
#   # cutree(hclust(dist(t(right.gx))), h = 0)
#   # 
#   # 
#   # View(ext.gx[,adj.groups == 3])
# }))
# 
# 
# hapstat = cbind(hapstat, suspicious[ match(1:nrow(hapstat), suspicious[,1]),])
# hapstat$extend[is.na(hapstat$extend)]=F
# hapstat$whichHaplo = 1:nrow(hapstat)












getHaplotypeSimilarirty = function(hapmatrix1, hapmatrix2){
  
  hapsim = setNames(as.data.frame(t(apply(expand.grid(1:ncol(hapmatrix1), 1:ncol(hapmatrix2)), 1, function(pair){
    h1 = hapmatrix1[,pair[1]]
    h2 = hapmatrix2[,pair[2]]
    percentMatch = sum(h1==h2, na.rm=T)/sum(!is.na(h1) & !is.na(h2))
    c(pair,percentMatch)
  }))), c("h1", "h2", "percentMatch"))
  
  hapsim=hapsim[order(hapsim$h1, hapsim$h2),]
  
  hapsim
  
}



checkPhaseOverlapWindown = function(fgt1,fgt2){
  
  
  hapsim = getHaplotypeSimilarirty(hapmatrix1=fgt1, hapmatrix2=fgt2)
  
  
  phase = do.call(cbind,lapply(c(1,3,5), function(n){
    matchg1=hapsim[hapsim$h1==n & hapsim$h2 %in% c(n,n+1),]
    matchg2=hapsim[hapsim$h1==n+1 & hapsim$h2 %in% c(n,n+1),]
    
    # Genome 1 of ovfi corresponds to genome 1 of ovfj 
    # &  genome 2 of ovfi corresponds to genome 2 of ovfj
    cis = matchg1$percentMatch[1] == 1 & matchg2$percentMatch[2] == 1
    
    # Genome 1 of ovfi corresponds to genome 2 of ovfj 
    # &  genome 2 of ovfi corresponds to genome 1 of ovfj
    trans = matchg1$percentMatch[2] == 1 & matchg2$percentMatch[1]  == 1
    
    c(cis=cis, trans=trans)
    
  }))
  
  colnames(phase) = data.table::tstrsplit(colnames(fi)[c(1,3,5)], ".g")[[1]]
  
  phase
  
}



findDivergentSnps = function(hapmatrix1, hapmatrix2){
  
  phasex=checkPhaseOverlapWindown(fgt1=hapmatrix1,fgt2=hapmatrix2)
  iscompatible = unlist(apply(phasex, 2, function(x){ x[1]|x[2]}))
  
  if(!sum(iscompatible)==length(iscompatible)){
    
    whichincomp = which(iscompatible==F)
    
    divergentsnps = unlist(lapply(whichincomp, function(wx){
      f1x = hapmatrix1[,(2*wx)-c(1,0)]
      f2x = hapmatrix2[,(2*wx)-c(1,0)]
      
      hsim=getHaplotypeSimilarirty(hapmatrix1=f1x, hapmatrix2=f2x)
      swap = sum(hsim$percentMatch[c(1,4)]) < sum(hsim$percentMatch[c(2,3)])
      if(swap){f2x = f2x[,c(2,1)]}
      divsnp = which(apply(f1x-f2x, 1, function(g){sum(g)!=0}))
      #divsnp = overlapwin[divsnp]
      divsnp
    }))
    
    divergentsnps=sort(unique(divergentsnps))
    return(divergentsnps)
  }else{
    return(NULL)
  }
  
}



InferFoundersHaploBlocks = function(founderA, founderB, founderM, rilsA, rilsB, info,nsnpWin = 100, maxSizeCM=3){
  
  
  
  windows = get.win(nrow(info), nsnpWin, nsnpWin)
  
  windows = do.call(c, lapply(1:length(windows), function(i){
    win = windows[[i]]
    gsize = diff(range(info$cM[win]))
    
    if (gsize < maxSizeCM) {
      # If within limits, return the first marker index of the window
      return(list(win))
    } else {
      print(win[i])
      # If not, return splitted window
      wheretosplit = which.min(abs(info$cM[win] - median(info$cM[win]) ))
      return(list(win[1:wheretosplit], win[(wheretosplit+1):length(win)]))
    }
    
  }))
  
  # Add intermediate windows which overlap with the adjacent windows
  windows = c(do.call(c, lapply(1:(length(windows)-1), function(i){
    wini = windows[[i]]
    winj = windows[[i+1]]
    midwin = c(wini[floor(length(wini)/2):length(wini)], winj[1:floor(length(winj)/2)])
    
    list(wini,midwin)
  })), list(windows[[length(windows)]]))
  
  
  
  founderhaplotypes = lapply(1:length(windows), function(step){
    #print(step)
    if(step %in% 10 == 0) print(step)
    
    # The window under consideration (this step to the next one -1)
    win =windows[[step]]
    
    # Infer founder haplotype blocks in the window + some stat about inference
    InferFoundersHaploBlocksWIN(win=win,
                                founderA=founderA, founderB=founderB, founderM=founderM,
                                rilsA=rilsA, rilsB=rilsB, 
                                info)
  })
  
  
  
  compatible = as.data.frame(do.call(rbind, lapply(1:(length(founderhaplotypes)-1), function(i){
    
    j=i+1
    fi = founderhaplotypes[[i]]$foundergt
    fj = founderhaplotypes[[j]]$foundergt
    wini = founderhaplotypes[[i]]$win
    winj = founderhaplotypes[[j]]$win
    
    winoverlap = wini[wini %in% winj]
    
    ovfi = fi[wini %in% winoverlap,]
    ovfj = fj[winj %in% winoverlap,]
    
    iscompatible = checkPhaseOverlapWindown(fgt1=ovfi, fgt2=ovfj)
    iscompatible = unlist(apply(iscompatible, 2, function(x){ x[1]|x[2]}))
    iscompatible
    
  })))
  
  incompatibility = which(apply(compatible,1,function(x){sum(x)<length(x)}) > 0)
  
  if( length(incompatibility)>0 ){
    print(incompatibility)
    
    if(sum(diff(incompatibility) <= 2)>0) stop("Code is probably not robust for the encoutered case")
    
    
    compatibleHaplo = lapply(incompatibility, function(i){
      wini=founderhaplotypes[[i]]$win
      winj=founderhaplotypes[[i+1]]$win
      winall = c(wini, winj)
      winall = sort(unique(winall))
      
      newfhap = InferFoundersHaploBlocksWIN(win=winall,
                                            founderA=founderA, founderB=founderB, founderM=founderM,
                                            rilsA=rilsA, rilsB=rilsB, 
                                            info)
      
      
      # Previously, block i was compatible with block i-1 & block j with j+1
      # So we want to check that we did not change that by with the reestimation of the founder haplotypes
      
      winh=founderhaplotypes[[i-1]]$win
      wink=founderhaplotypes[[i+2]]$win
      
      
      iscompatible = cbind(#Check left compatibility
        checkPhaseOverlapWindown(fgt1=founderhaplotypes[[i-1]]$foundergt[winh %in% newfhap$win,], 
                                 fgt2=newfhap$foundergt[newfhap$win %in%  winh,]),
        #Check right compatibility
        checkPhaseOverlapWindown(fgt1=founderhaplotypes[[i+2]]$foundergt[wink %in% newfhap$win,], 
                                 fgt2=newfhap$foundergt[newfhap$win %in%  wink,]))
      
      iscompatible = c(unlist(apply(iscompatible, 2, function(x){ x[1]|x[2]})))
      iscompatible = sum(iscompatible)==length(iscompatible)
      
      if(!iscompatible) stop("Code needs to be improved here")
      
      return(newfhap)
      
    })
    
    
    for(i in 1:length(incompatibility)){
      founderhaplotypes[[ incompatibility[i] ]] <- compatibleHaplo[[i]]
    }
    
    founderhaplotypes = founderhaplotypes[-(incompatibility+1)]
    
  }
  
  ## Get rid of the overlaping windows
  #founderhaplotypes = founderhaplotypes[c(1:(length(founderhaplotypes)-1) %% 2 == 1, T)]
  
  return(founderhaplotypes)
}



highdistwin=which(unlist(lapply(founderhaplotypes, function(x){ x$totdist/length(x$win) }))>0)
splithere=c(0,which(diff(highdistwin)>1))
highdistwin=lapply(1:(length(splithere)-1), function(i){
  highdistwin[(splithere[i]+1):(splithere[i+1])]
})


i=75

whichwin = highdistwin[[i]]
extwhichwin=c(whichwin[1]-1, whichwin, whichwin[length(whichwin)]+1)
#extwhichwin = whichwin

highdisthap = founderhaplotypes[ extwhichwin ]
extwin=sort(unique(unlist(lapply(highdisthap, function(x){x$win}))))
#extwin = extwin[-which(extwin==33850)]
sum(unique(unlist(lapply(highdisthap, function(x){x$totdist}))))
sum(unique(unlist(lapply(highdisthap, function(x){x$nbreaks}))))

exthap = InferFoundersHaploBlocksWIN(win=extwin,
                                     founderA=founderA, founderB=founderB, founderM=founderM,
                                     rilsA=rilsA, rilsB=rilsB, 
                                     info)


diff(DivergentsSnps)

DivergentsSnps = sort(unique(unlist(lapply(highdisthap, function(x){
  
  overlapwin = exthap$win[exthap$win %in% x$win]
  f1 = exthap$foundergt[exthap$win %in% overlapwin,]
  f2 = x$foundergt[x$win %in% overlapwin,]
  
  divsnps = findDivergentSnps(hapmatrix1=f1, hapmatrix2=f2)
  overlapwin[divsnps]
  
}))))


if(length(DivergentsSnps)>0){
  extwin = extwin[-which(extwin %in% DivergentsSnps)]
  exthap = InferFoundersHaploBlocksWIN(win=extwin,
                                       founderA=founderA, founderB=founderB, founderM=founderM,
                                       rilsA=rilsA, rilsB=rilsB, 
                                       info)
  DivergentsSnps = sort(unlist(lapply(highdisthap, function(x){
    
    overlapwin = exthap$win[exthap$win %in% x$win]
    f1 = exthap$foundergt[exthap$win %in% overlapwin,]
    f2 = x$foundergt[x$win %in% overlapwin,]
    
    divsnps = findDivergentSnps(hapmatrix1=f1, hapmatrix2=f2)
    overlapwin[divsnps]
    
  })))
  
}




i=675
j=705

fi = founderhaplotypes[[i]]$foundergt
wini = founderhaplotypes[[i]]$win
fj = founderhaplotypes[[j]]$foundergt
winj = founderhaplotypes[[j]]$win

gi = cbind(fi, cbind(rilsA, rilsB)[wini,])
gj = cbind(fj, cbind(rilsA, rilsB)[winj,])


groupsi = cutree(hclust(dist(t(gi))), h = 0) 
table(groupsi)

groupsj = cutree(hclust(dist(t(gj))), h = 0) 
table(groupsj)

fcomb = table(groupsi[groupsi %in% groupsi[1:5] & groupsj %in% groupsj[1:4]], 
              groupsj[groupsi %in% groupsi[1:5] & groupsj %in% groupsj[1:4]])

c()

View(cbind(groupsi, 
           groupsj))

swap=which.max(c(sum(fcomb[cbind(c(1,2), c(1,2))]),  #cis numbers
                 sum(fcomb[cbind(c(1,2), c(2,1))]))) # trans numbers






((f1[,3]+f1[,4])/2)-((founderB[overlapwin,1]+founderB[overlapwin,2])/2)
((f2[,3]+f2[,4])/2)-((founderB[overlapwin,1]+founderB[overlapwin,2])/2)


p=f1[,4]
pj=f2[,4]

gx = cbind(f1[,4],f2[,4], rilsB[overlapwin,])


groups = cutree(hclust(dist(t(gx))), h = 0) 
table(groups)
fcomb = table(groupsi[groupsi %in% groupsi[1:2] & groupsj %in% groupsj[1:2]], 
              groupsj[groupsi %in% groupsi[1:2] & groupsj %in% groupsj[1:2]])

swap=which.max(c(sum(fcomb[cbind(c(1,2), c(1,2))]),  #cis numbers
                 sum(fcomb[cbind(c(1,2), c(2,1))]))) # trans numbers



#heterozygosity = cbind(whichwin=1:length(founderhaplotypes), do.call(rbind, lapply(founderhaplotypes, function(x){x$heterozygous})))


SWAP = as.data.frame(do.call(rbind, lapply(1:(length(founderhaplotypes)-1), function(i){
  #print(i)
  j=i+1
  fi = founderhaplotypes[[i]]$foundergt
  fj = founderhaplotypes[[j]]$foundergt
  
  wini = founderhaplotypes[[i]]$win
  winj = founderhaplotypes[[j]]$win
  
  winoverlap = wini[wini %in% winj]
  
  ovfi = fi[wini %in% winoverlap,]
  ovfj = fj[winj %in% winoverlap,]
  
  phasepair = checkPhaseOverlapWindown(fgt1=ovfi, fgt2=ovfj)
  #phasepair look like this:
  #       FA    FB   FM
  #cis    TRUE  TRUE TRUE
  #trans FALSE FALSE TRUE
  # cis=T: the two windows are already in the same phase
  # trans=T: the two windows are in opposite phase
  # cis & trans = T; happen when fi & fj are homozygous within 
  # cis & trans = F, means that code is not robust to some case and need to be revised
  
  # Let's transform the cis/trans matrix in a "swap vector" where T = swap (trans); F = no swap (cis)
  # 
  sapply(1:ncol(phasepair), function(n){
    swapx = c(F,T)[which(phasepair[,n])]
    
    # If the two swap are possible, check if fj is homozygous
    # if fj is hom, the phase does not matter so swap = T or F would be good, we choose F
    # In other case (ex: fi/fj are hom within winoverlap but fj is still heterozygous winthin winh ):
    # => we cannot infer the phase like this
    # => And will need to use other info later
    if(length(swapx)==2 & founderhaplotypes[[j]]$het[n]){swapx=NA}
    if(length(swapx)==2 & !founderhaplotypes[[j]]$het[n]){swapx=F}
    if(length(swapx)==0) stop("This should not happen. Need some revision")
    
    swapx
  })
  
})))

SWAP = rbind(rep(F,3), SWAP)

#29401,29450
#589, 590
which(lapply(founderhaplotypes, function(x){x$win[1]}) == 29450)

heterozygosity = do.call(rbind, lapply(founderhaplotypes, function(x){x$heterozygous}))

phasedfounderhaplotypes = lapply(1:ncol(SWAP), function(thisfounder){
  
  swapx = SWAP[,thisfounder]
  breaks = c(1,which(is.na(swapx)), length(swapx))
  
  PHASED = do.call(rbind,lapply(1:(length(breaks)-1), function(b){
    print(b)
    wins = breaks[b]:(breaks[b+1]-1)
    subswapx = swapx[wins]
    phasex = rep(1, length(wins))
    
    for(swaphere in which(subswapx)){
      phasex[swaphere:length(phasex)] = phasex[swaphere:length(phasex)]*-1
    }
    
    ##Just ensure that each "phased block" start with phase 1, event though, I think it is already the case
    ## Should be ok even if it is not 
    #phasex=phasex*phasex[1]
    
    phasedhaplo = lapply(1:length(wins), function(whichwin){
      x = founderhaplotypes[[ wins[whichwin] ]]
      whichgt = (thisfounder*2)-c(1,0)
      if(phasex[whichwin] == -1){ whichgt=whichgt[c(2,1)] }
      list(win=x$win,foundergt=x$foundergt[,whichgt], gsize=x$gsize)
    })
    
    #Now we want to fuse haplotypes knowing there is an overlap (we want to delete)
    # Will will start deleting overlaping from the second window of this block. We store the first window here.
    firsthap = cbind(block=b,whichsnp=phasedhaplo[[1]]$win, phasedhaplo[[1]]$foundergt) 
    
    if(length(phasedhaplo)>2) wid = 2:length(phasedhaplo)
    if(length(phasedhaplo)==2) wid = 2
    if(length(phasedhaplo)<2) wid = NULL
    
    
    if(!is.null(wid)){
      phasedhaplo=do.call(rbind, lapply(wid, function(whichwin){
        print(whichwin)
        pi= phasedhaplo[[ whichwin-1 ]]$foundergt
        pj = phasedhaplo[[ whichwin ]]$foundergt
        wini = phasedhaplo[[ whichwin-1 ]]$win
        winj = phasedhaplo[[ whichwin ]]$win
        
        winoverlap = wini[wini %in% winj]
        
        #Verification that the window are in good phase to be sure (their overlaping genotype should be the same)
        if(sum(abs(c(pi[wini %in% winoverlap,] - pj[winj %in% winoverlap,]))) > 0){
          stop("Windows not correctly phased. Should not happens. There is a bug.")}
        
        
        cbind(block=b,whichsnp=winj[!(winj %in% winoverlap)], pj[!(winj %in% winoverlap),])
      }))
      
      
      phasedhaplo=rbind(firsthap,phasedhaplo)
      
    }else{
      
      phasedhaplo=firsthap
      
    }
    
    
    
  }))
  
})


names(phasedfounderhaplotypes) = c("FA", "FB", "FM")

thisfounder = 2

phap = as.data.frame(phasedfounderhaplotypes[[1]])

i=33
j=i+1
pi=phap[phap$block==i,]
pj=phap[phap$block==j,]

if(names(phasedfounderhaplotypes)[thisfounder]=="FA") rx = rilsA
if(names(phasedfounderhaplotypes)[thisfounder]=="FB") rx = rilsB
if(names(phasedfounderhaplotypes)[thisfounder]=="FM") rx = cbind(rilsA,rilsB)

gi = cbind(pi[,3:4], rx[pi$whichsnp,])
gj = cbind(pj[,3:4], rx[pj$whichsnp,])

groupsi = cutree(hclust(dist(t(gi))), h = 0) 
groupsj = cutree(hclust(dist(t(gj))), h = 0) 


fcomb = table(groupsi[groupsi %in% groupsi[1:2] & groupsj %in% groupsj[1:2]], 
              groupsj[groupsi %in% groupsi[1:2] & groupsj %in% groupsj[1:2]])

swap=which.max(c(sum(fcomb[cbind(c(1,2), c(1,2))]),  #cis numbers
                 sum(fcomb[cbind(c(1,2), c(2,1))]))) # trans numbers


cbind(pi[,3:4], founderB[pi$whichsnp,])


swap=ifelse(swap==1,F,T)
swap



InferFoundersHaploBlocksWIN(win=pi$whichsnp,
                            founderA=founderA, founderB=founderB, founderM=founderM,
                            rilsA=rilsA, rilsB=rilsB, 
                            info)


