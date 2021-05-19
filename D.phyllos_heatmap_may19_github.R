##Phyllosphere genetic comparison
#Acer VanWallendael

library(readr)
library(tidyverse)
library(vegan)
library(readxl)

setwd("~/Desktop/phyllos_all/rdata")
#Start with distance matrix fungi @ kbs
load("~/Desktop/phyllos_all/rdata/otus_dist_bray2_jan27.rda")
#metadata
load("~/Desktop/phyllos_all/rdata/OTU_funkb_jan27.rda")
load("../rdata/OTU_fungi_jan27.rda")

ps.funkbs<-OTU_funkb
meta_samp <- as.data.frame(sample_data(ps.funkbs))
colnames(meta_samp)<-c("PLOT_GL", "Date","Site", "inf","ID","Contam")
meta_samp$Code<-(ps.funkbs@sam_data@row.names)

#Get genetic distance. Generated with data from Lovell et al. 2021 Nature
load("~/Desktop/GWAS/rdata/gwas_mat_all.rda")
gwas_mat_all$X786<-NULL
gwas_all<-colnames(gwas_mat_all)

#give plant_id to otutab
dist_mat1<-as.matrix(otus_dist_bray2)
DOE <- read_excel("~/Downloads/DOE_GWAS_Master Plant List (2).xlsx")
plantids<-DOE$PLANT_ID[match(substr(rownames(dist_mat1),1,5), DOE$PLOT_GL)]
rownames(dist_mat1)==meta_samp$Code

#split by time
plant_idlist<-split(plantids, ps.funkbs@sam_data$Date)

shared<-intersect(plant_idlist$`158`,intersect(plant_idlist$`212`,intersect(plant_idlist$`233`, intersect(plant_idlist$`260`, intersect(gwas_all,plant_idlist$`286`)))))
#check if all shared samples are present
length(shared)==92

#sub the genetic dist mat
gwas_mat_sub<-gwas_mat_all[gwas_all %in% shared,
                           gwas_all %in% shared]

genmat2<-as.matrix(gwas_mat_sub)

heatmap(genmat2)

#need to set the order
distance    = dist(genmat2)
cluster     = hclust(distance, method="ward.D")
dendrogram  = as.dendrogram(cluster)
Rowv        = rowMeans(genmat2, na.rm = T)
dendrogram  = reorder(dendrogram, Rowv)
reorderfun = function(d,w) { d }
heatmap(genmat2,Colv=dendrogram ,Rowv=dendrogram, reorderfun=reorderfun)

otu_allib<-ps.funkbs@otu_table@.Data
meta_samp$PLANT_ID<-DOE$PLANT_ID[match(meta_samp$PLOT_GL, DOE$PLOT_GL)]
otus_names<-meta_samp$PLANT_ID

calc_heatmap<-function(otutab, date, meta_samp, otus_names, genmat2){
  otu_allib<-otutab
#sep allib by time
time1<-otu_allib[,meta_samp$Date==date]
time1_plant<-otus_names[meta_samp$Date==date]
time1_plant2<-time1_plant[!duplicated(time1_plant)]
#rm duplicates
time1_<-time1[,!duplicated(time1_plant)]
time1_order<-as.data.frame(time1_[,time1_plant2 %in% shared])
ncol(time1_order)==length(shared)

otus_rel1<-decostand(t(time1_order), method="hellinger", na.rm=T)
otus_dist1<-vegdist(otus_rel1, method = "bray", na.rm = T)
otus_distmat1<-as.matrix(otus_dist1, labels=T)

#check order
colnames(otus_distmat1)<-meta_samp$PLANT_ID[match(colnames(otus_distmat1),meta_samp$Code)]
rownames(otus_distmat1)<-meta_samp$PLANT_ID[match(rownames(otus_distmat1),meta_samp$Code)]
otus_dist_ord1<-otus_distmat1[colnames(genmat2),colnames(genmat2)]
return(otus_dist_ord1)
}

dates<-list(158,212,233,260,286)

dists<-lapply(dates, FUN = function(x){calc_heatmap(otu_allib,x,meta_samp, otus_names,genmat2)})

#PLOT####

#pdf(file = "../figs/heatmaps_jan27.pdf", height=5, width=5)
heatmap(genmat2,Colv=dendrogram ,Rowv=dendrogram, reorderfun=reorderfun, main = "Genetic Distance")
heatmap(1-dists[[1]],Colv=dendrogram ,Rowv=dendrogram, reorderfun=reorderfun, main = "Community Distance DOY 158")
heatmap(1-dists[[2]],Colv=dendrogram ,Rowv=dendrogram, reorderfun=reorderfun, main = "Community Distance DOY 212")
heatmap(1-dists[[3]],Colv=dendrogram ,Rowv=dendrogram, reorderfun=reorderfun, main = "Community Distance DOY 233")
heatmap(1-dists[[4]],Colv=dendrogram ,Rowv=dendrogram, reorderfun=reorderfun, main = "Community Distance DOY 260")
heatmap(1-dists[[5]],Colv=dendrogram ,Rowv=dendrogram, reorderfun=reorderfun, main = "Community Distance DOY 286")
#graphics.off()

names(dists)<-c("Time1_158", "Time2_212", "Time3_233", "Time4_260", "Time5_286")
save(dists, file="dists_heatmap_jan27.rda")

t1<-mantel((1-dists[[1]]), genmat2)
t2<-mantel((1-dists[[2]]), genmat2)
t3<-mantel((1-dists[[3]]), genmat2)
t4<-mantel((1-dists[[4]]), genmat2)
t5<-mantel((1-dists[[5]]), genmat2)

mantels<-data.frame(rval=c(t1$statistic,t2$statistic,t3$statistic,t4$statistic,t5$statistic),
                    DOY=c(158,212,233,260,286))

B<-ggplot(mantels, aes(x=DOY, y=rval))+
  geom_point(col="navyblue")+
  geom_line(col="navyblue")+
  geom_text(label=round(mantels$rval,3),nudge_y = c(.04,-.02,.02,.015,.025),
            nudge_x = c(3,5,-5,0,-3))+
  labs(x="Day of Year (DOY)", 
       # title="Genetic-Microbial community correlations",
       y="Mantel's R")+
  theme_classic()
B
ggsave("../figs/mantel_heatmap_jan27.png", height = 4, width=4.5)
ggsave("../figs/mantel_heatmap_jan27.pdf", height = 4, width=4.5)
ggsave("../figs/mantel_heatmap_jan27.eps", height = 4, width=4.5)
save.image("~/Desktop/phyllos_all/rdata/D.phyllos_heatmap_jan27.RData")

#fig for paper####
library(patchwork)
heatmap(genmat2,Colv=dendrogram ,Rowv=dendrogram, reorderfun=reorderfun, main = "Genetic Distance")
heatmap(1-dists[[1]],Colv=dendrogram ,Rowv=dendrogram, reorderfun=reorderfun, main = "Community Distance DOY 158")
heatmap(1-dists[[4]],Colv=dendrogram ,Rowv=dendrogram, reorderfun=reorderfun, main = "Community Distance DOY 260")
B

