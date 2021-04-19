library(tidyverse)
library(bedr)
library(biomaRt)
library(tidyr)
library(jsonlite)


args = commandArgs(trailingOnly=TRUE)

#load data
height <- read_tsv(args[1])#"50_irnt.gwas.imputed_v3.both_sexes.tsv"
blocks <- read_tsv(args[2],col_names=c("chr","start","end"))#"fourier_ls-chr3.bed"
variants <- read_tsv(args[3])#"variants.tsv"

#subset chr 3 
height_chr3 <- filter(height, str_detect(variant, "^3:"))
variants_chr3 <- filter(variants, str_detect(variant, "^3:"))

# split the "variant" column into separate columns
height_chr3 <- tidyr::separate(height_chr3, "variant", c("chr", "pos", "major", "minor"))
height_chr3$pos<-as.numeric(height_chr3$pos)

#subset only to those columns of interest
variants_chr3<-variants_chr3 %>% dplyr::select(pos,ref,alt,rsid,varid)

# make a new column that has the compact notation for the bed regions
height_chr3 <- mutate(height_chr3, bed=paste0("chr", chr,":", pos-1, "-", pos))
blocks <- mutate(blocks, bed=paste0(chr,":", start, "-", end))

# join the LD blocks onto the list of SNPs
  regions <- bedr.join.region(height_chr3$bed, blocks$bed) %>%
  mutate(ld_block=paste0(V4, ":", V5, "-", V6)) %>%
    select(index, ld_block)

# now attach this to the original df
height_chr3 <- inner_join(height_chr3, regions, by = c("bed"="index"))

# get the best p-value per ld_block 
independentA <- height_chr3 %>% 
  group_by(ld_block) %>% 
  slice_min(pval, n = 1)

independentB <- height_chr3 %>% 
  filter(pval <= 5e-8) %>%
  group_by(ld_block) %>% 
  slice_min(pval, n = 1)

# subset only those variants that have rsID
independent_shortA<-inner_join(variants_chr3, independentA, by = c("pos"="pos"))
independent_shortB<-inner_join(variants_chr3, independentB, by = c("pos"="pos"))

# get columns of interest from the variants file  and filter those that have proper rsID 
independent_shortA<-inner_join(variants_chr3, independentA, by = c("pos"="pos")) %>%
                    filter(str_detect(rsid, "^rs"))
independent_shortB<-inner_join(variants_chr3, independentB, by = c("pos"="pos")) %>%
                    filter(str_detect(rsid, "^rs"))

# retrieve ancestral allele
database_url<-"http://grch37.rest.ensembl.org/variation/human/"
ancestral_alleleA<-c()
summary_data<-data.frame()
for( i in 1:length(independent_shortA$rsid)){
  message("Retrieving id \t", i)
  summary_data<-fromJSON(paste0(database_url,independent_shortA$rsid[i]))
  ancestral_alleleA[i]<-summary_data$mappings$ancestral_allele
}

# polarize by ancestral allele 
for(i in 1:length(independent_shortA$rsid)){
  if(is.na(ancestral_alleleA[i])){
    ancestral_alleleA[i]<-independent_shortA$major[i]
  }
  if(isTRUE(independent_shortA$major[i] != ancestral_alleleA[i])){ # change sign 
    independent_shortA$beta[i]<-(-1)*independent_shortA$beta[i]
  }
}
