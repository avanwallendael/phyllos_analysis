##Phyllosphere analyses - NMDS & PERMANOVA
#Acer VanWallendael

library(readxl)
library(tidyverse)
library(vegan)
library(phyloseq)
library(statip)

setwd("~/Desktop/phyllos_all/rdata")

#Data Setup#####
#load decontaminated scaled fungi phyloseq obj
#and the kbs-specific obj
 load("../rdata/OTU_fungi_jan27.rda")
 load("../rdata/OTU_funkb_jan27.rda")
load('../rdata/OTU_shared_jan27.rda')

#KBS MDS analysis####
ps.samp<-OTU_funkb

#check
unique(ps.samp@sam_data@.Data[[5]])

#get metadata from ps
meta_samp <- as.data.frame(do.call(cbind,ps.samp@sam_data@.Data))
colnames(meta_samp)<-c("PLOT_GL", "Date","Site", "inf","ID","PLANT_ID", "Eco","Subpopulation", "Contam")
meta_samp$Code<-(ps.samp@sam_data@row.names)

#split out otus
otus<-ps.samp@otu_table@.Data

otus_rel<-decostand(t(otus), method="hellinger")
otus_dist_bray2<-vegdist(otus_rel, method = "bray")
betadisper1 <- betadisper(as.dist(otus_rel), group=meta_samp$Code)
nmds2<-metaMDS(otus_rel, k=3, trymax=100)

plot(nmds2, "sites")
orditorp(nmds2, "sites")
#orditorp(nmds1, "species")
plot(betadisper1)
pcoa_beta<-betadisper1$vectors

#prep for plot
points2<-nmds2$points
points2<-as.data.frame(points2)
points2$Code<-meta_samp$Code
points3<-left_join(points2, meta_samp, by="Code")

nmds2$stress
#0.1024153

#Site differences####
#split out otus
otus_site<-OTU_shared@otu_table@.Data

meta_samp2 <- as.data.frame(do.call(cbind,OTU_shared@sam_data@.Data))
colnames(meta_samp2)<-c("PLOT_GL", "Date","Site", "inf","ID","PLANT_ID", "Eco","Subpopulation", "Contam", "idnew")
meta_samp2$Code<-(OTU_shared@sam_data@row.names)

otus_site_rel<-decostand(t(otus_site), method="hellinger")
otus_dist_bray3<-vegdist(otus_site_rel, method = "bray")
nmds_site<-metaMDS(otus_site_rel, k=3, trymax = 100)

plot(nmds_site, "sites")
orditorp(nmds_site, "sites")
#orditorp(nmds_site, "species")

#prep for plot
points_site<-nmds_site$points
points_site<-as.data.frame(points_site)
points_site$Code<-meta_samp2$Code
points_site2<-left_join(points_site, meta_samp2, by="Code")

nmds_site$stress
#0.1034034

#Distances and Permanova####
#can't have NAs in phenos. 
meta_samp$Eco[is.na(meta_samp$ECO)]<-"Intermediate"
meta_samp$numdate<-as.numeric(meta_samp$Date)

#full model KBS
full_margin<-adonis2(as.dist(otus_dist_bray2)~
        numdate+Subpopulation+inf, 
        data=meta_samp,
        by="margin")
full_margin

#try with day1 as term
meta_samp$day1<-points3$MDS2[match(meta_samp$PLOT_GL, points3$PLOT_GL)]

full_margin2<-adonis2(as.dist(otus_dist_bray2)~
                       day1+numdate+inf+Subpopulation, 
                     data=meta_samp,
                     by="margin")
full_margin2
#need to use adonis() if you have interaction term
int_terms<-adonis(as.dist(otus_dist_bray2)~
                      numdate:Subpopulation+inf, 
        data=meta_samp)
int_terms

full_terms<-adonis2(as.dist(otus_dist_bray2)~
                       numdate+Subpopulation+numdate:Subpopulation+inf, 
                     data=meta_samp,
                     by="terms")
full_terms

#full model all
meta_samp2$numdate<-as.numeric(meta_samp2$Date)

all_margin<-adonis2(as.dist(otus_dist_bray3)~
                       Site+numdate+inf+Subpopulation, 
                     data=meta_samp2,
                     by="margin")
all_margin
all_terms<-adonis2(as.dist(otus_dist_bray3)~
                     Site+numdate+Subpopulation+inf, 
                    data=meta_samp2,
                    by="terms")
all_terms

adonis_tests<- list(full_margin, full_terms, all_margin, all_terms)

# put into a better output
adonis_tests2<-vector("list", 4)
for(i in 1:4){
  adonis_tests2[[i]]<-as.data.frame(adonis_tests[[i]])
  #adonis_tests2[[i]][nrow(adonis_tests2[[i]]),ncol(adonis_tests2[[i]])]<-attr(adonis_tests2[[i]], "heading")[2]
  adonis_tests2[[i]]$SumOfSqs<-round(adonis_tests2[[i]]$SumOfSqs,1)
  adonis_tests2[[i]]$R2<-round(adonis_tests2[[i]]$R2,3)
  adonis_tests2[[i]]$`F`<-round(adonis_tests2[[i]]$`F`,2)
  rownames(adonis_tests2[[i]])[rownames(adonis_tests2[[i]])=="numdate"]<-"DOY"
  rownames(adonis_tests2[[i]])[rownames(adonis_tests2[[i]])=="numdate:Subpopulation"]<-"DOY:Subpopulation"
  rownames(adonis_tests2[[i]])[rownames(adonis_tests2[[i]])=="inf"]<-"Infection"
  }

library(kableExtra)
options(knitr.kable.NA = '')

adonis_tests2[[1]] %>%
  kbl(caption = "KBS Margin") %>%
  kable_classic(full_width = F, html_font = "Cambria")%>%
  row_spec((1:(nrow(adonis_tests2[[1]])-2)),bold=T,hline_after = F)%>%
save_kable("../figs/Permanova_table_KBS_margin_feb1.pdf")

adonis_tests2[[2]] %>%
  kbl(caption = "KBS Sequential") %>%
  kable_classic(full_width = F, html_font = "Cambria")%>%
  row_spec((1:(nrow(adonis_tests2[[2]])-2)),bold=T,hline_after = F)%>%
  save_kable("../figs/Permanova_table_KBS_terms_mar17.pdf")

adonis_tests2[[3]] %>%
  kbl(caption = "Multi-site Margin") %>%
  kable_classic(full_width = F, html_font = "Cambria")%>%
  row_spec((1:(nrow(adonis_tests2[[3]])-2)),bold=T,hline_after = F)%>%
  save_kable("../figs/Permanova_table_multi_margin_feb1.pdf")

adonis_tests2[[4]] %>%
  kbl(caption = "Multi-Site Sequential") %>%
  kable_classic(full_width = F, html_font = "Cambria")%>%
  row_spec((1:(nrow(adonis_tests2[[4]])-2)),bold=T,hline_after = F)%>%
  save_kable("../figs/Permanova_table_multi_terms_feb1.pdf")

#alltests<-do.call(rbind, adonis_tests2)
#alltests[6,5]
#save output
save(nmds2, file=                               "nmds2_jan27.rda")
save(nmds2, file=                           "nmds_site_jan27.rda")
save(points3, file=                      "points3_otus_jan27.rda")
save(points_site2, file =                "points_site2_jan27.rda")
write_csv(alltests, file =                   "alltests_jan27.csv")
save(otus_dist_bray2, file=  "../rdata/otus_dist_bray2_jan27.rda")
save(otus_dist_bray3, file=  "../rdata/otus_dist_bray3_jan27.rda")
save.image("~/Desktop/phyllos_all/rdata/B.phyllos_nmds_jan27.RData")