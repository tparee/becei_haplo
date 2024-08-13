
#####################################
#### UPLOAD FUNCTION AND DATA  ######

# Function:
source("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/haplotype/RILs_haplotypes/function_findhaplo_in_RILS.R")

# Founder phased haplotypes:
load(file="/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/haplotype/founderHaplotypes/240515_beceiFounderPhasedHaplotypes_chrV.Rdata")

# RILs genotype Tables
rils <- read_csv("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/genotype/StringentFiltering/becei_chr_V_Cross_AB_Founder_RIL_GT_table.csv")


######################
## Prepare data ######
mm = match(phasedfounderhaplotypes$POS, rils$POS) # Match positions 
rils=rils[mm,]

snps = phasedfounderhaplotypes[,1:4] # snps info
colnames(snps)=tolower(colnames(snps))

founders = phasedfounderhaplotypes[,5:10] # Keep only genotype info
foudernames = colnames(founders) # extract foundernames
founders = matrix(as.numeric(unlist(c(founders ))), ncol=ncol(founders)) # format as numeric matrix
colnames(founders)=foudernames

# Split founder of cross A and founder of cross B
foundersA = founders[,grepl("A|M", colnames(founders))]
foundersB = founders[,grepl("B|M", colnames(founders))]

# Format rils as numric matrix
rils = rils[,grepl("_QG", colnames(rils))]
rilnames =  colnames(rils)
rils[rils=="1/1"] = "1"
rils[rils=="0/0"] = "0"
rils[rils=="0/1"]=NA
rils = matrix(as.numeric(unlist(c(rils))), ncol=ncol(rils))
colnames(rils) = rilnames

# Split RILs of cross A and founder of cross B
rilsA = rils[, grepl("A_Q",colnames(rils))]
rilsB = rils[, grepl("B_Q",colnames(rils))]



# Infer haplotype blocks for rils A
ix = 1:ncol(rilsA)
RIlsA_FounderHaplotypeBlocks_chrV = do.call(rbind,parallel::mclapply(ix,mc.cores = 2, function(thisril){
  print(thisril)
  ril = rilsA[,thisril]
  ril[ril!=1 & ril!=0]=NA
  output = try(haplosearch(ril, founders=foundersA, snps))
  if(class(output)!='try-error'){
    output$rilname = colnames(rilsA)[thisril]
    return(output)
  }else{
    return(NULL)
  }
}))




#RIlsB_FounderHaplotypeBlocks_chrIII$pos1=snps$pos[RIlsB_FounderHaplotypeBlocks_chrIII$whichsnp1]
#RIlsB_FounderHaplotypeBlocks_chrIII$pos2=snps$pos[RIlsB_FounderHaplotypeBlocks_chrIII$whichsnp2]
#RIlsB_FounderHaplotypeBlocks_chrIII$cM1=snps$cm[RIlsB_FounderHaplotypeBlocks_chrIII$whichsnp1]
#RIlsB_FounderHaplotypeBlocks_chrIII$cM2=snps$cm[RIlsB_FounderHaplotypeBlocks_chrIII$whichsnp2]


#RIlsB_FounderHaplotypeBlocks$rilname = colnames(rilsB)[match(RIlsB_FounderHaplotypeBlocks$rilname, colnames(rils))]

#save(RIlsB_FounderHaplotypeBlocks, file="~/Documents/Documents - MacBook Pro de tom/rockmanlab/cbecei_haplo/240504_rilsBfounderHaplotypeBlocks_chrI.Rdata")


#save(RIlsB_FounderHaplotypeBlocks_chrIII, file="~/Documents/Documents - MacBook Pro de tom/rockmanlab/cbecei_haplo/240505_rilsBfounderHaplotypeBlocks_chrIII.Rdata")


# Infer haplotype blocks for rils B
ix = 1:ncol(rilsB)
RIlsB_FounderHaplotypeBlocks_chrV = do.call(rbind,parallel::mclapply(ix,mc.cores = 2, function(thisril){
  print(thisril)
  ril = rilsB[,thisril]
  ril[ril!=1 & ril!=0]=NA
  output = try(haplosearch(ril, founders=foundersB, snps))
  if(class(output)!='try-error'){
    output$rilname = colnames(rilsB)[thisril]
    return(output)
  }else{
    return(NULL)
  }
}))



RIlsA_FounderHaplotypeBlocks_chrV$pos1=snps$pos[RIlsA_FounderHaplotypeBlocks_chrV$whichsnp1]
RIlsA_FounderHaplotypeBlocks_chrV$pos2=snps$pos[RIlsA_FounderHaplotypeBlocks_chrV$whichsnp2]
RIlsA_FounderHaplotypeBlocks_chrV$cM1=snps$cm[RIlsA_FounderHaplotypeBlocks_chrV$whichsnp1]
RIlsA_FounderHaplotypeBlocks_chrV$cM2=snps$cm[RIlsA_FounderHaplotypeBlocks_chrV$whichsnp2]

save(RIlsA_FounderHaplotypeBlocks_chrV, file="~/Documents/Documents - MacBook Pro de tom/rockmanlab/cbecei_haplo/240516_rilsAfounderHaplotypeBlocks_chrV.Rdata")


RIlsB_FounderHaplotypeBlocks_chrV$pos1=snps$pos[RIlsB_FounderHaplotypeBlocks_chrV$whichsnp1]
RIlsB_FounderHaplotypeBlocks_chrV$pos2=snps$pos[RIlsB_FounderHaplotypeBlocks_chrV$whichsnp2]
RIlsB_FounderHaplotypeBlocks_chrV$cM1=snps$cm[RIlsB_FounderHaplotypeBlocks_chrV$whichsnp1]
RIlsB_FounderHaplotypeBlocks_chrV$cM2=snps$cm[RIlsB_FounderHaplotypeBlocks_chrV$whichsnp2]

save(RIlsB_FounderHaplotypeBlocks_chrV, file="~/Documents/Documents - MacBook Pro de tom/rockmanlab/cbecei_haplo/240516_rilsBfounderHaplotypeBlocks_chrV.Rdata")
