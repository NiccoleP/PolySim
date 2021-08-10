library(tidyr)
library(tidyverse)
library(cowplot)

vcf_list <- list.files(pattern = "*.vcf")
###Load all files
for(i in 1:length(vcf_list)){
  filepath <- file.path("/home/niccole/Documents/R_POLY_SIM/",paste(vcf_list[i]))
  assign(vcf_list[i], read.delim(filepath,sep = "\t"))
}

generations<-str_extract(vcf_list,"[0-9]+(?=\\.)")

pos_eff<-read.table("/home/niccole/position_effect.txt")
colnames(pos_eff)<-c("POS","effect")
pos_eff$POS<-pos_eff$POS+1

vcf_files<-mget(ls(pattern=".vcf"))

for (i in 1:length(vcf_files)) {
    vcf_files[[i]]<-merge(vcf_files[[i]],pos_eff,by.x='POS',by.y='POS',all.x = TRUE)
}

vcf_copy<-vcf_files

for (i in 1:length(vcf_files)) {
  vcf_files[[i]]<-vcf_files[[i]][, 10:length(vcf_files[[i]])]
}

for (i in 1:(length(vcf_files))){
  vcf_files[[i]][,]<-case_when(vcf_files[[i]][,] == '0|0' ~ 0,
                     vcf_files[[i]][,] == '0|1' ~ 1,
                     vcf_files[[i]][,] == '1|0' ~ 1,
                     vcf_files[[i]][,] == '1|1' ~ 2)
}


for (i in 1:length(vcf_files)) {
  vcf_files[[i]]$effect<-vcf_copy[[i]]$effect
}


for (i in 1:length(vcf_files)) {
vcf_files[[i]][,]<-vcf_files[[i]][,]*vcf_files[[i]]$effect
}

for (i in 1:length(vcf_files)) {
  vcf_files[[i]]$effect<-NULL
}

vcf_means<-list()
for (i in 1:length(vcf_files)) {
  vcf_means[[i]]<-vcf_files[[i]] %>% summarise(across(starts_with("i"), mean))
}

for (i in 1:length(vcf_means)) {
  vcf_means[[i]]<-data.frame(unlist(vcf_means[[i]]))
  }

for (i in 1:length(vcf_means)) {
 colnames(vcf_means[[i]])<-c("means")
}


for (i in 1:length(vcf_means)) {
  vcf_means[[i]]$generation<-rep(generations[i],length(vcf_means[[1]]$means))
}

all_p0<-bind_rows(vcf_means)

all_p0$generation<-as.numeric(all_p0$generation)

sorted_p0<-all_p0 %>% arrange(all_p0$generation)

gen_mean_p0<-tibble(aggregate(sorted_p0[,1], list(sorted_p0$generation), mean))
colnames(gen_mean_p0)<-c("gen","mean")

cols<-c("African pop (YRI)"="indianred1","European pop (CEU)"="goldenrod1","Asian pop (CHB)"="lightblue4")

ggplot(sorted_p0,aes(x=generation,y=means))+xlab("generation")+ylab("polygenic score")+
  geom_jitter(alpha=0.1,color='indianred')+
  geom_jitter(data=sorted_p1,aes(x=generation,y=means),alpha=0.5,color='goldenrod1')+ ## I repeated this same script for another population and store the df as sorted_p1 and sorted_p2
  geom_jitter(data=sorted_p2,aes(x=generation,y=means),alpha=0.5,color='lightblue4')+
  geom_line(data=gen_mean_p0,aes(x=gen,y=mean), color="firebrick4")+
  geom_line(data=gen_mean_p1,aes(x=gen,y=mean), color="darkorange")+
  geom_line(data=gen_mean_p2,aes(x=gen,y=mean) ,color="gray24")+
  geom_hline(yintercept = -0.03782096,color="firebrick4")+
  geom_hline(yintercept = -0.0541799,color="darkorange")+
  geom_hline(yintercept = -0.05647234,color="gray24")+
  theme_cowplot()
