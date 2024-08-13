######################
### HAPLO pipeline ###

(!) WIP. For now it is still very dirty and not optimal.
(!) The genotype folder is too big to be added on Git Hub. It should contain two folder: "stringentFiltering" and "Unfiltered", containing the high quality SNP and RILs genotype table, or unfiltered ones, respectively. 

This pipeline aims to infer the founder haplotypes and filter the genotyping error in becei RILS

=> STEP 1: infer the founder haplotypes (AUTOSOMES ONLY FOR NOW): 

Use script: beceiFounders_phasing.R

In short, we start from the "genotype/StringentFiltering/becei_chr_#_Cross_AB_Founder_RIL_GT_table.csv".
This files contains genotype table of high quality SNV for filtered RILs and poolseq inferred genotypes.
The code look at the segregating haplotypes in RILs to infer the founder haplotypes.
The founder haplotypes are saved as "beceiFounderPhasedHaplotypes_chr#.Rdata"

Done for: I, II, III, IV, V
Missing: X (code not ready)



=> STEP 2: Infer the founder haplotype blocks in (good quality) rils: 

Use script: RILs_haplotypes/searchFounderHaploInBeceiRils.R

In this step, we find the founder haplotype block within the RILs, for each cross (A or B) separately.
The file is saved as "yymmdd_rils#founderHaplotypeBlocks_chr#.Rdata".

the "240502_plot_becei.R" can be used to visualize haplotype in RILs


Done for: I, II, III, IV, V
Missing: V, X




=> STEP 3: Add the filtered out SNVs

Use scripts: getFormat2.R & incorporate_missing_snps.R

In this step, we use the yymmdd_rils#founderHaplotypeBlocks_chr#.Rdata generated in STEP 2 to test each SNV that was previously filtered out.
First, we need a little change of format which is done with script getFormat2.R. and saved as rilsfounderHaplotypeBlocks_chr#_format2.Rdata.
Then, the script incorporate_missing_snps.R use the latter file to icorporate filtered out SNV.
This is done at each missing SNV that they match with the haplotypes blocks, i.e., they are monomorphic within each haplotype block origin. 
The filtered out SNV we try to incorporate are found in genotype/Unfiltered/chr#_genotypeTable_beceiRILs_unfiltered.csv
(!) this dataset was filtered

Done for: I (!) the dataset was not totally unfiltered, so it should be redone with a truly unfiltered dataset
Missing:  II, III, IV, V, X

=> STEP 4: Check wether the heterozygosity make sense given the founder haplotype

This is done by trying to match a given RIL genotype with the a combination of founder haplotypes (hom: founderA.g1 + founderA.g1;  or het: founderA.g1 + founderM.g1) in 1,000 SNV window.
Only the minimal genotype distance is kept. If heterozygosity is legit, the (minimal) genotype distance should remain low (basal level).

For chromosome I, most heterozygosity does not seems legit. Maybe a few case it is. 

This step may be by-passed for now, and just consider every heterozygous SNV is NA.



=> STEP 4 BIS: re-infering founder haplotype blocks but with all the new SNV we incorporated

For now I would consider every single SNV as NA.

and then it is re-runing STEP 2, but with the new dataset and for all the RILs.

Then, knowing the founder haplotypes, it should be relatively easy to write a small piece of code that infer the the missing SNPs.

ALTERNATIVELY, this STEP could be replace by other software that are designed to infer missing genotypes with known founder haplotypes. 





SUMMARY TABLE:

               I     II    III    IV    V   X
STEP 1         x      x     x      x    x
STEP 2         x      x     x      x
STEP 3        (x)     
STEP 4
STEP 4 BIS








 
