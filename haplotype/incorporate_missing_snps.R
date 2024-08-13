library(data.table)
library(readr)

fNAmax=0.1 # max proportion of NA value for a snp
genodistMax = 2 # max genotype distance from the expected genotype given the assumed founder haplotype


source("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/haplotype/RILs_haplotypes/function_findhaplo_in_RILS.R")
load("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/haplotype/founderHaplotypes/240509_beceiFounderPhasedHaplotypes_chrI.Rdata")
load("~/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/haplotype/RILs_haplotypes/rilsfounderHaplotypeBlocks_chrI_format2.Rdata")



# extract info of the snps of the phased founder 
# i.e., the stringently filtered snps
goodsnps = phasedfounderhaplotypes[,1:4]
goodsnps$pos=goodsnps$POS

# Now let's import the RILs genotype table with all the snps
rils_unfiltered <- read_csv("~/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/genotype/Unfiltered/chrI_genotypeTable_beceiRILs_unfiltered.csv")
snps_unfiltered = rils_unfiltered[,1:6]

# Just keep all the lines in which we identify the founder haplotype blocks
rils_unfiltered=rils_unfiltered[,colnames(rils_unfiltered) %in%  c(rilsHaplo$rilname,rilsHaplo$rilname)]

# identify the snps that were previously filtered out
# i.e., snps that are in snps_unfiltered but absent in goodsnps
nomatch = which(is.na(match(snps_unfiltered$id, goodsnps$ID)))




supplFounderGeno = do.call(rbind, lapply(nomatch, function(i){
  
  # for each snps that were previously filtered out
  # we are going to see if there are consistent with the known linkage blocks
  
  if(i %% 1000 == 0) print(i)
  
  ppos = snps_unfiltered[i,]$pos
  
  startpos = min(rilsHaplo$pos1)
  endpos = max(rilsHaplo$pos2)
  
  if(ppos<startpos){ppos = startpos}
  if(ppos>endpos){ppos = endpos}
  
  # All the haplotype blocks among rils that encompass the target position
  blocks = rilsHaplo[(rilsHaplo$pos1 <= ppos & rilsHaplo$pos2 >= ppos),]
  
  # Add the genotype observed in the rils at the target position
  blocks$genotype = unlist(c(rils_unfiltered[i,match(blocks$rilname, colnames(rils_unfiltered))]))
  
  # Verify there is not too much NA
  if( sum(is.na(blocks$genotype))/nrow(blocks) > fNAmax ) return(NULL)
  
  # Verify that there is only one allele per founder haplotype 
  nperhap = table(blocks$founder, blocks$genotype)
  founders = strsplit(rownames(nperhap),";")
  genovalue = as.numeric(colnames(nperhap))
  
  # Ensure there is variation 0/1 at this snp (if only 0/0.5, not good)
  if(sum(c(0,1) %in% genovalue)<2) return(NULL) 
  isuniquefounder = unlist(lapply(founders, length))==1
  
  # A
  foundersgeno = do.call(rbind, lapply(which(isuniquefounder), function(f){
    
    x=nperhap[f,] # genotype distribution for this founder
    wgeno = which.max(x) # The most frequent genotype = the genotype of this founder
    gdist=sum((genovalue*x)[-wgeno]) # number of individual with a different genotype than expected
    
    data.frame(founders=rownames(nperhap)[f], genotype = genovalue[wgeno], genodist=gdist)
    
  }))
  
  
  # Sum the genodist observed in haplotype blocks corresponding to a unique founder
  # and the ones from the blocks corresponding to blocks which can correspond do several founders
  genotypedistance = sum(foundersgeno$genodist) + sum(unlist(lapply(which(!isuniquefounder), function(f){
    
    x=nperhap[f,]
    found = founders[[f]]
    
    isbad = !(names(x) %in% as.character(foundersgeno[match(found, foundersgeno$founders),]$genotype))
    (genovalue*x)[isbad]
    
  })))
  
  # if genotype distance over threshold, snp is bad
  if(genotypedistance > genodistMax){ return(NULL)}
  
  # if snp is good, return the ganotype for each founder
  foundersgeno = foundersgeno$genotype[match(c("FA.g1","FA.g2","FB.g1","FB.g2","FM.g1","FM.g2"),foundersgeno$founders)]
  if(sum(is.na(foundersgeno))>0){return(NULL)}
  if(sum(foundersgeno) %in% c(0,length(foundersgeno)) | sum(foundersgeno==0.5)>0){return(NULL)}
  
  foundersgeno = c(i, foundersgeno)
  
  
}))


#x=suppFounderGeno[1,2:7]

#keep = apply(suppFounderGeno[,2:7],1, function(x){
#  
#  if(sum(is.na(x))>0){return(F)}
#  if( sum(x) %in% c(0,length(x)) | sum(x==0.5)>0){return(F)}else{return(T)} 
#  
#  })

#suppFounderGeno=suppFounderGeno[keep,]
supplFounderGeno=as.data.frame(supplFounderGeno)
colnames(supplFounderGeno) = c("whichsnp", "FA.g1","FA.g2","FB.g1","FB.g2","FM.g1","FM.g2")
supplFounderGeno = cbind(snps_unfiltered[supplFounderGeno$whichsnp,], supplFounderGeno[,2:7])


mm=match(tolower(colnames(phasedfounderhaplotypes)), tolower(colnames(supplFounderGeno)))
mm2 = mm
mm2[is.na(mm)]=1
supplFounderGeno=supplFounderGeno[,mm2]
supplFounderGeno[,is.na(mm)]=NA
colnames(supplFounderGeno)=colnames(phasedfounderhaplotypes)

phasedfounderhaplotypes = rbind(phasedfounderhaplotypes, supplFounderGeno)
phasedfounderhaplotypes=phasedfounderhaplotypes[,-which(colnames(phasedfounderhaplotypes)=="genotypedistance")]
phasedfounderhaplotypes$cM = approx(x=phasedfounderhaplotypes$POS, y=phasedfounderhaplotypes$cM, xout = phasedfounderhaplotypes$POS)$y

phasedfounderhaplotypes=phasedfounderhaplotypes[order(phasedfounderhaplotypes$POS),]
save(phasedfounderhaplotypes, file="/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/haplotype/founderHaplotypes/beceiFounderPhasedHaplotypes_chrI_complete.Rdata")






