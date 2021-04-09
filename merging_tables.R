library(httr)
library(jsonlite)
library(xml2)
library(ggplot2)
library(MASS)
library(fitdistrplus)
library(tidyverse)
library("bedr")

LD_C3_AFR<-read.table("/home/niccole/Documents/RACIMOLAB_RESOURCES/FILES/AFR_ls_chr3.bed")
C3_50_irnt<-read.table("/home/niccole/Documents/RACIMOLAB_RESOURCES/FILES/gwas_50_irnt_chrom3_pos_major")
UKBB_3<-read.table("/home/niccole/Documents/RACIMOLAB_RESOURCES/FILES/UKBB_3.info",header = TRUE)
#UKBB_3[2]<-NULL
#C3_50_irnt[1]<-NULL
colnames(LD_C3_AFR)<-c("chr","start","end")
LD_C3_AFR$ld_segment<-c(1:185)
colnames(C3_50_irnt)<-c("minor_allele","minor_AF","low_confidence_variant","n_complete_samples","AC","ytx","beta","se","tstat","pval","POS","major_allele")
C3_50_irnt<-C3_50_irnt[c(11,12,1,2,7,8,10)]
C3_50_irnt_RS_merged<-inner_join(UKBB_3,C3_50_irnt,by=c("POS"))

gwas_50_irnt_c3.bed<-data.frame(rep(paste0("chr3"),948269),C3_50_irnt_RS_merged$POS,C3_50_irnt_RS_merged$POS,C3_50_irnt_RS_merged$RS)
colnames(gwas_50_irnt_c3.bed)<-c("chr","start","end","rs")

a.sort.n<-bedr.sort.region(gwas_50_irnt_c3.bed, method = "natural")
b.sort.n<-bedr.sort.region(LD_C3_AFR, method = "natural")
num_intersec_LD_50_irnt<-bedr.join.region(a.sort.n,b.sort.n)
summary_LD<-data.frame(num_intersec_LD_50_irnt$start,num_intersec_LD_50_irnt$rs,num_intersec_LD_50_irnt$ld_segment)
colnames(summary_LD)<-c("POS","id","ld_segment")
summary_50_irnt<-C3_50_irnt_RS_merged[c(1,3,4,5,6,9,10)]#rs pos major minor MAF MAF2 beta

AFR_LD_and_qtl_info<-inner_join(summary_50_irnt,summary_LD,by=c("POS"))
FINAL_TABLE<-AFR_LD_and_qtl_info %>% mutate(Diff = POS -lag(POS))
rs<- 'rs'
FINAL<-subset(FINAL_TABLE,grepl(rs,RS))
