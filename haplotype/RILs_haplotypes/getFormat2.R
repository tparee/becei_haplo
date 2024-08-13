library(data.table)
library(readr)
source("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/haplotype/RILs_haplotypes/function_findhaplo_in_RILS.R")
load("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/haplotype/founderHaplotypes/240509_beceiFounderPhasedHaplotypes_chrI.Rdata")
load("~/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/haplotype/RILs_haplotypes/240504_rilsAfounderHaplotypeBlocks_chrI.Rdata")
load("~/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/haplotype/RILs_haplotypes/240504_rilsBfounderHaplotypeBlocks_chrI.Rdata")



#RIlsA_FounderHaplotypeBlocks$whichsnp1 = sapply(RIlsA_FounderHaplotypeBlocks$pos1, function(x){max(which(x>=phasedfounderhaplotypes$POS))})
#RIlsA_FounderHaplotypeBlocks$whichsnp2 = sapply(RIlsA_FounderHaplotypeBlocks$pos2, function(x){min(which(x<=phasedfounderhaplotypes$POS))})
#RIlsB_FounderHaplotypeBlocks$whichsnp1 = sapply(RIlsB_FounderHaplotypeBlocks$pos1, function(x){max(which(x>=phasedfounderhaplotypes$POS))})
#RIlsB_FounderHaplotypeBlocks$whichsnp2 = sapply(RIlsB_FounderHaplotypeBlocks$pos2, function(x){min(which(x<=phasedfounderhaplotypes$POS))})


snps = phasedfounderhaplotypes[,1:4]
snps$pos=snps$POS

RIlsA_FounderHaplotypeBlocks=do.call(rbind,lapply(split(RIlsA_FounderHaplotypeBlocks, RIlsA_FounderHaplotypeBlocks$rilname), function(HAPLO){
  
  #HAPLO=split(RIlsA_FounderHaplotypeBlocks, RIlsA_FounderHaplotypeBlocks$rilname)[[3]]
  
  print(HAPLO$rilname[1])
  HAPLO = FuseHapSeparatedByOneSNP(HAPLO=HAPLO, posinfo = c("whichsnp","pos", "cM"),
                                   foundercolnames=c("FA.g1", "FA.g2", "FM.g1", "FM.g2"))
  
  
  HAPLO=divideOverlapHapBlocksByFounderID(HAPLO, foundercolnames=c("FA.g1", "FA.g2", "FM.g1", "FM.g2"),
                                          info=snps)
  
  
  HAPLO
  
}))


RIlsB_FounderHaplotypeBlocks=do.call(rbind,lapply(split(RIlsB_FounderHaplotypeBlocks, RIlsB_FounderHaplotypeBlocks$rilname), function(HAPLO){
  
  
  print(HAPLO$rilname[1])
  HAPLO = FuseHapSeparatedByOneSNP(HAPLO=HAPLO, posinfo = c("whichsnp","pos", "cM"),
                                   foundercolnames=c("FB.g1", "FB.g2", "FM.g1", "FM.g2"))
  
  
  HAPLO=divideOverlapHapBlocksByFounderID(HAPLO, foundercolnames=c("FB.g1", "FB.g2", "FM.g1", "FM.g2"),
                                          info=snps)
  
  
  HAPLO
  
}))



rilsHaplo = rbind(RIlsA_FounderHaplotypeBlocks, RIlsB_FounderHaplotypeBlocks)

save(rilsHaplo,file = "~/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/haplotype/RILs_haplotypes/rilsfounderHaplotypeBlocks_chrI_format2.Rdata")
