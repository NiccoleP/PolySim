library(vcfR)
library(tidyr)
library(tidyverse)
library(jsonlite)
library(argparser)

# read data 
p <- arg_parser("optimal phenotype calculation ")
p <- add_argument(p, "pop_1", help="population one")
p <- add_argument(p, "pop_2", help="population two")
p <- add_argument(p, "pop_3", help="population three")
p <- add_argument(p, "var_nopol", help="independent variants before polarization")
args <- parse_args(p)

# load data

pop_1<- read.table(args$pop_1)
#pop_1_indivs_and_variants.dataframe
pop_1<-pop_1%>%filter(str_detect(V5, "^rs"))
pop_2<- read.table(args$pop_2)
#pop_2_indivs_and_variants.dataframe
pop_2<-pop_2%>%filter(str_detect(V5, "^rs"))
pop_3<- read.table(args$pop_3)
#pop_3_indivs_and_variants.dataframe
pop_3<-pop_3%>%filter(str_detect(V5, "^rs"))
var_nopol<-read.table(args$var_nopol)
#variants_before_polarizing
colnames(var_nopol)<-c("pos","ref","alt","rsid","beta","ancestral")
all_pops<-list(pop_1,pop_2,pop_3)

for( i in 1:length(all_pops)){
  all_pops[[i]]$beta<-var_nopol$beta
}

# retrieve ancestral allele
database_url<-"http://grch37.rest.ensembl.org/variation/human/"
ancestral_allele<-c()
summary_data<-data.frame()
for( i in 1:length(all_pops[[1]]$V5)){
  message("Retrieving id \t", i)
  summary_data<-fromJSON(paste0(database_url,all_pops[[1]]$V5[i]))
  ancestral_allele[i]<-summary_data$mappings$ancestral_allele
}

#Ancestral allele must be present in the biallelic calls 
for(i in 1:length(all_pops)){
all_pops[[i]]<-subset(all_pops[[i]],all_pops[[i]]$V3==ancestral_allele | all_pops[[i]]$V4 == ancestral_allele )
}


# polarize by ancestral allele 
## CHB 
for(j in 1:length(all_pops)){
for(i in 1:length(all_pops[[1]]$V5)){
  if(is.na(ancestral_allele[i])){
    ancestral_allele[i]<-all_pops[[j]]$V3[i] 
  }
  if(isTRUE(all_pops[[j]]$V3[i] != ancestral_allele[i])){ # change sign 
    all_pops[[j]]$beta[i]<-(-1)*all_pops[[j]]$beta[i]
  }
}
}

tmp<-matrix(0,ncol = ncol(all_pops[[1]]), nrow = nrow(all_pops[[1]]))
tmp2<-matrix(0,ncol = ncol(all_pops[[2]]), nrow = nrow(all_pops[[2]]))
tmp3<-matrix(0,ncol = ncol(all_pops[[3]]), nrow = nrow(all_pops[[3]]))

tmps<-list(tmp,tmp2,tmp3)

for (k in 1:length(all_pops)) {
for(j in 6:ncol(all_pops[[k]])){
  for(i in 1:length(all_pops[[k]]$V1)){
    if(isFALSE(all_pops[[k]]$V4[i] == ancestral_allele[i])){ # alternative is derived 
      message("alt is not ancestral", i)
      if(isTRUE(all_pops[[k]][i,j]=='0|0')){
        message("ALT == derived then 0|0 is 0 ", i )
        tmps[[k]][i,j]<-+0
      }
      if(isTRUE(all_pops[[k]][i,j]=='0|1'|all_pops[[k]][i,j]=='1|0' )){
        message("ALT == derived then 1|0 or 0|1 is 1 ", i )
        tmps[[k]][i,j]<-+1
      }
      if(isTRUE(all_pops[[k]][i,j]=='1|1')){
        message("ALT == derived then 1|1 is 2 ", i )
        tmps[[k]][i,j]<-+2
      }
    }
    
    else{message("alt is ancestral", i)  # alternative is not derived 
      if(isTRUE(all_pops[[k]][i,j]=='0|0')){
        message("ref == derived then 0|0 is 2", i )
        tmps[[k]][i,j]<-+2
      }
      if(isTRUE(all_pops[[k]][i,j]=='0|1'|all_pops[[k]][i,j]=='1|0' )){
        message("ref == derived then 1|0 or 0|1 is 1 ", i )
        tmps[[k]][i,j]<-+1
      }
      if(isTRUE(all_pops[[k]][i,j]=='1|1')){
        message("ref == derived then 1|1 is 0 ", i )
        tmps[[k]][i,j]<-+0
      }
    }
  }
}
}


for(i in 1:length(tmps)){
  tmps[[i]]<-tmps[[i]][,-c(1:5,length(all_pops[[i]]))]
}

for(i in 1:length(tmps)){
  tmps[[i]]<-tmps[[i]]*all_pops[[i]]$beta
}

summary_stats<-c(c("pop_1",mean(colSums(tmps[[1]])),sd(colSums(tmps[[1]]))),
c("pop_2",mean(colSums(tmps[[2]])),sd(colSums(tmps[[2]]))),
c("pop_3",mean(colSums(tmps[[3]])),sd(colSums(tmps[[3]]))))
write_delim(x =as.data.frame(summary_stats),file = "summary_stats_pheno",delim = '\t',col_names = FALSE)
