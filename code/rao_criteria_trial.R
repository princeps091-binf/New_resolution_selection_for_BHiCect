library(tidyverse)
library(Matrix)
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


gg_res_scatter<-chr_dat_l[[5]] %>% ggplot(.,aes(KR,raw))+geom_point()+scale_x_log10()+scale_y_log10()

ggsave("./img/KR_vs_raw_10kb.png",gg_res_scatter)
do.call(bind_rows,chr_dat_l[1:4])%>% 
  ggplot(.,aes(KR,raw))+
  geom_point()+
  facet_wrap(res~.)+
  scale_x_log10()+scale_y_log10()
chr_dat_l[[3]] %>% mutate(IO=raw>=1000) %>% summarise(Rao.criteria=sum(IO)/n())
