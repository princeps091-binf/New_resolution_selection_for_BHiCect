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
neighbour_HiC_tbl<-do.call(bind_rows,lapply(chr_dat_l,function(dat){
  
  dat %>% filter(X1!=X2) %>% mutate(d=abs(X1-X2)) %>% filter(d<=min(d))
})) %>% 
  mutate(res=fct_relevel(res,res_set)) 

neighbour_HiC_tbl%>% 
  ggplot(.,aes(raw,color=res))+geom_density()+
  facet_wrap(res~.,scales="free")

neighbour_HiC_tbl %>% 
  group_by(res) %>% 
  summarise(m=sum(raw>=1000)/n())

neighbour_HiC_tbl %>% 
  group_by(res) %>% 
  summarise(m=median(raw))
