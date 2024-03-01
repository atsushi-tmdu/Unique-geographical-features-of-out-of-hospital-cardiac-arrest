library(readxl)
library(dplyr)
library(ggplot2)
library(maptools)
library(INLA)
library(data.table)
library(GISTools)
library(tidyverse)
library(sf)
library(pROC)
sf_use_s2()

df<- read_xlsx("../pre_analysys/all7.xlsx")
jpnmap2<-rgdal::readOGR("./preanalysis/limited_map_without_lesions/tokyomap.shp")
df$longitude= df$longtitude
legionmap <-rgdal::readOGR("./preanalysis/limited_map_with_lesions/tokyomap.shp")
inla_df <- df %>% dplyr::select(id,gender,year,age,byCPR,vf_vt,cpc,time1,time2,time3,time4,witness,latitude,longitude) %>% filter(vf_vt==1)
inla_df <- inla_df %>%  mutate(cpc1_2 =if_else(cpc<3, true=1, false=0))
inla_df <- inla_df %>% mutate(delta_time1=time2-time1, delta_time2=time3-time2,delta_time3=time4-time3)
inla_df <- inla_df %>% mutate(id=row_number())
inla_df<-inla_df %>% dplyr::select(id,gender,year,age,byCPR,cpc1_2,delta_time1,delta_time2,delta_time3,witness,,latitude,longitude)
inla_df<-mutate(inla_df, total_time = if_else(condition = delta_time1 > 0 & delta_time2 > 0 & delta_time3 > 0,
                                    true = delta_time1 + delta_time2 + delta_time3,
                                    false = -1))
inla_df<-inla_df %>% filter(total_time!=-1)
inla_df$normalized_total_time<- as.numeric(inla_df$total_time)/as.numeric(inla_df$total_time %>% mean())

#geographical
coords<- inla_df %>%  dplyr::select(longitude,latitude)
coords$admYM<-NULL
pj<-CRS("+proj=longlat +ellps=WGS84 +no_defs")
jpnmap2<-spTransform(jpnmap2,CRS=pj)
colnames(coords) <- c("x","y")
coordinates(coords) <- ~x+y
sp <-SpatialPoints(coords, proj4string = pj)
legionmap<-spTransform(legionmap,CRS=pj)


dt2_list_df<-inla_df %>% dplyr::select(delta_time2,cpc1_2) %>% arrange(delta_time2) %>% mutate(cum_cpc1_2=cumsum(cpc1_2)) %>% group_by(delta_time2) %>% summarise(cum_cpc1_2=max(cum_cpc1_2),n=n()) %>% mutate(cum_n=cumsum(n)) %>% ungroup() %>% mutate(ratio=cum_cpc1_2/cum_n)
dt2_list_df %>% dplyr::filter(delta_time2<13*60) %>% mutate(delta_time2_min=(delta_time2/60))
model<-lm(data=dt2_list_df %>% dplyr::filter(delta_time2<13*60),ratio ~ delta_time2)



g1<-ggplot(data=dt2_list_df,aes(x=delta_time2/60,y=ratio*100,weight=cum_n))+
  geom_smooth()+
  theme_bw()+
  ylim(20,80)+
  xlab('Time spent on site (min)')+
  ylab('Patients with CPC 1-2 (%)')
ggsave(file="./time_and_mortality/delta_time2.pdf",g1, dpi=800)


dt3_list_df<-inla_df %>% dplyr::select(delta_time3,cpc1_2) %>% arrange(delta_time3) %>% mutate(cum_cpc1_2=cumsum(cpc1_2)) %>% group_by(delta_time3) %>% summarise(cum_cpc1_2=max(cum_cpc1_2),n=n()) %>% mutate(cum_n=cumsum(n)) %>% ungroup() %>% mutate(ratio=cum_cpc1_2/cum_n)

g2<-ggplot(data=dt3_list_df,aes(x=delta_time3,y=ratio*100,weight=cum_n*100))+
  geom_smooth()+
  theme_bw()+
  ylim(20,80)+
  xlab('Time spent for transfer (min)')+
  ylab('Patients with CPC 1-2  (%)')
ggsave(file="./time_and_mortality/delta_time3.pdf",g2, dpi=800)