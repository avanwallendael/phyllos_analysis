##Phyllosphere figures
#Acer VanWallendael

setwd("~/Desktop/phyllos_all/figs")

library(tidyverse)
library(vegan)
library(readxl)

#NMDS plot####
load("~/Desktop/phyllos_all/rdata/points3_otus_jan27.rda")
load("~/Desktop/phyllos_all/rdata/nmds2_jan27.rda")
#DOE <- read_excel("~/Downloads/DOE_GWAS_Master Plant List (2).xlsx")
PHEN_META_TRIM <- read_csv("~/Downloads/PHEN_META_TRIM (3).csv")

load("~/Desktop/phyllos_all/rdata/OTU_funkb_jan27.rda")
ps.samp<-OTU_funkb

points1<-points3
points1$Lat<-PHEN_META_TRIM$LATITUDE[match(points1$PLOT_GL, PHEN_META_TRIM$PLOT_GL)]
points1$Lon<-PHEN_META_TRIM$LONGITUDE[match(points1$PLOT_GL, PHEN_META_TRIM$PLOT_GL)]

#write for gwas
write_csv(points1, "../rdata/points1_kbs_fungi_jan6.csv")

stress<-round(nmds2$stress,3)

A<-ggplot(points1)+
  geom_point(aes(x=MDS1, y=MDS2, fill=Date, shape=Subpopulation),col="black", size=3, na.rm=T)+
  scale_fill_manual( values=c("#8435a0","#c2a5cf","white", "#a6dba0","#008837"))+
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  #scale_size_manual(values=c(1,2,3,4,5))+
  scale_shape_manual(values = c(21,22,  23,24,23), name="Subpopulation")+
  #facet_grid(.~Date, space = "free", scales = "free")+
  annotate("text", x = -1.2,y=-0.45, label=paste("Stress = ", stress), size=6)+ 
  labs(x="NMDS1",y="NMDS2", title="a")+
   theme_classic(base_size = 12)
A
ggsave("../figs/Pop_otu_phyllos_feb28.png", height=4, width = 9)
ggsave("../figs/Pop_otu_phyllos_feb28.pdf", height=4, width = 9)
ggsave("../figs/Pop_otu_phyllos_feb28.eps", height=4, width = 9)


#Site differences####
load("../rdata/points_site2_jan27.rda")
load("../rdata/OTU_fungi_jan27.rda")

phenologykey<-data.frame(date=unique(points_site2$Date),
                         season=c("Mid", "Mid", "Early", "Mid", "Mid", "Late","Late", "Late", "Early","Early", "Mid", "Late","Late", "Late", "Early"))
points_site2$phenology<-phenologykey$season[match(points_site2$Date, phenologykey$date)]
points_site2$Site<-factor(points_site2$Site, levels = c("M", "C", "P","K"))
points_site2$Site<-c( "KBS", "Columbia", "Austin", "Kingsville")[match(points_site2$Site, c( "M", "C", "P","K"))]
points_site2$Site<-factor(points_site2$Site, levels = c( "KBS", "Columbia", "Austin", "Kingsville"))
points_site2$Phenology<-factor(points_site2$phenology, levels = c("Early", "Mid", "Late"))

site<-ggplot(points_site2)+
  geom_point(aes(x=MDS1, y=MDS2, col=Phenology , shape=Site), na.rm=T, size=2, stroke=1)+
  #geom_line(aes(x=-MDS1, y=MDS2, group=PLOT_GL))+
  #facet_grid(Subpopulation~.)+
  scale_color_manual(values=c( "#7b3294","#d5b8dc","#008837"))+
  annotate("text", x=-.7, y=.6, label="Northern", size=4.5)+
  annotate("text", x=1, y=-.55, label="Southern", size=4.5)+
  labs(y="NMDS2", x="NMDS1")+
  scale_shape_manual(values = c(4,8,0,2), name="Site")+
  theme_classic(base_size = 12)

site
ggsave("../figs/nmds_site_feb28.png", height=3, width=5)
ggsave("../figs/nmds_site_feb28.pdf", height=3, width=5)
ggsave("../figs/nmds_site_feb28.eps", height=3, width=5)


#Map####
library(maps)

map1<-map_data("usa")
states_map <- map_data("state")

KBSM <- PHEN_META_TRIM
kbsm2<-KBSM[which(KBSM$LONGITUDE<(-40)),]

table(kbsm2$SUBPOP_SNP)

subpops<-unique(kbsm2$SUBPOP_SNP)
bigpops<-c("Atlantic", "Atlantic", "Gulf", "Midwest", "Gulf", "Gulf", "Gulf", "Undetermined", "Midwest")
kbsm2$Subpopulations<-bigpops[match(kbsm2$SUBPOP_SNP, subpops)]

kbsm3 <- kbsm2 %>% filter(!is.na(Subpopulations))%>%filter(!is.na(ECOTYPE_SNP_CHLR))

ggplot(map1) +
  geom_polygon(aes(long, lat, group = group), color = "black", fill="white")+
  geom_point(data=kbsm3, aes(x=LONGITUDE, y=LATITUDE, col=Subpopulations, shape=ECOTYPE_SNP_CHLR), size=1.5)+
  labs(x="Longitude", y="Latitude", title="Original Collection Locations")+
  scale_color_manual(name="Genetic Subpopulation",values = c("dodgerblue", "salmon", "darkorchid4"))+
  scale_shape(name="Morphological Ecotype")+
  theme_classic()
#ggsave("map_phyllospheresites_nov12.png", height=3, width = 6)
#ggsave("map_phyllospheresites_nov12.pdf", height=3, width = 6)
#ggsave("map_phyllospheresites_nov12.eps", height=3, width = 6)


