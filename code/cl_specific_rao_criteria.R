# Examine whether the clusters found with previous criteria fit the Rao criteria
library(tidyverse)
library(data.tree)
options(scipen=999999999)
#############################################################
res_set<-c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num<-c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
KR_dat_file<-'~/Documents/multires_bhicect/data/GM12878/'
raw_dat_file<-'~/Documents/multires_bhicect/data/GM12878/raw/'

chromo<-"chr22"

chr_dat_l<-lapply(res_set,function(x){
  
  KR_dat<-read_delim(file = paste0(KR_dat_file,x,'/',chromo,'.txt'),delim = '\t',col_names = F)
  raw_dat<-read_delim(file = paste0(raw_dat_file,x,'/',chromo,'.txt'),delim = '\t',col_names = F)
  chr_dat<-KR_dat %>% dplyr::rename(KR=X3) %>% left_join(.,raw_dat%>% dplyr::rename(raw=X3)) %>% mutate(res=x)
  
})

chr_dat_l<-lapply(chr_dat_l,function(x){
  x%>%filter(!(is.nan(KR)))
})
names(chr_dat_l)<-res_set
#############################################################
## Load the original BHiCect results
spec_res_file<-'~/Documents/multires_bhicect/data/GM12878/spec_res/'
base::load(paste0(spec_res_file,chromo,'_spec_res.Rda'))
tmp_cl<-names(chr_spec_res$cl_member)[sample(1:length(chr_spec_res$cl_member),1)]
tmp_cl<-names(chr_spec_res$cl_member)[1]

cl_res<-strsplit(tmp_cl,split="_")[[1]][1]
cl_tbl<-tibble(chr=chromo,res=cl_res,cl=tmp_cl,bins=list(as.numeric(chr_spec_res$cl_member[[tmp_cl]]))) %>% 
  unnest(cols=c(bins))
chr_dat_l[[unique(cl_tbl$res)]] %>% 
  filter(X1 %in% cl_tbl$bins | X1 %in% cl_tbl$bins) %>% 
  summarise(Rao.criteria=quantile(raw,0.2)>=1000,quant.1k=sum(raw>=1000)/n())
#############################################################
## Load the original BHiCect results 
spec_res_file<-'~/Documents/multires_bhicect/data/GM12878/spec_res/Rao_criteria/'
base::load(paste0(spec_res_file,chromo,'_spec_res.Rda'))
tmp_cl<-names(chr_spec_res$cl_member)[sample(1:length(chr_spec_res$cl_member),1)]
#tmp_cl<-names(chr_spec_res$cl_member)[1]

cl_res<-strsplit(tmp_cl,split="_")[[1]][1]
cl_tbl<-tibble(chr=chromo,res=cl_res,cl=tmp_cl,bins=list(as.numeric(chr_spec_res$cl_member[[tmp_cl]]))) %>% 
  unnest(cols=c(bins))
chr_dat_l[[unique(cl_tbl$res)]] %>% 
  filter(X1 %in% cl_tbl$bins | X1 %in% cl_tbl$bins) %>% 
  summarise(Rao.criteria=quantile(raw,0.2),quant.1k=sum(raw>=1000)/n())
