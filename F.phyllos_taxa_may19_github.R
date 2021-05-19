#Testing models for ecological patterns.
#Acer VanWallendael

library(vegan)
library(tidyverse)

setwd("~/Desktop/phyllos_all")
#lead w OTU_funkb
load("~/Desktop/phyllos_all/rdata/OTU_funkb_jan31.rda")


#time series####
#Try from trajectories in vegclust https://cran.r-project.org/web/packages/vegclust/vignettes/CTA.html
library(vegclust)
library(RColorBrewer)
library(smacof)
library(MASS)


#first subsetted distmat.
load("~/Desktop/phyllos_all/rdata/otus_dist_bray2_jan27.rda")
meta<-as.data.frame(as.matrix(sample_data(OTU_funkb)))
#check that they line up
attr(otus_dist_bray2, "Labels")[1:5]==rownames(meta)[1:5]

#split by population
otu_mid<-prune_samples(meta$Subpopulations=="Midwest", OTU_funkb)
otu_gul<-prune_samples(meta$Subpopulations=="Gulf", OTU_funkb)
otu_atl<-prune_samples(meta$Subpopulations=="Atlantic", OTU_funkb)
otu_int<-prune_samples(meta$Subpopulations=="Intermediate", OTU_funkb)

#pull tabs & clean
middf <- as.data.frame(t(otu_table(otu_mid)))
guldf <- as.data.frame(t(otu_table(otu_gul)))
atldf <- as.data.frame(t(otu_table(otu_atl)))
intdf <- as.data.frame(t(otu_table(otu_int)))

middf$PLOT_GL<-substr(rownames(middf),1,5)
guldf$PLOT_GL<-substr(rownames(guldf),1,5)
atldf$PLOT_GL<-substr(rownames(atldf),1,5)
intdf$PLOT_GL<-substr(rownames(intdf),1,5)

middf$Date<-substr(rownames(middf),7,9)
guldf$Date<-substr(rownames(guldf),7,9)
atldf$Date<-substr(rownames(atldf),7,9)
intdf$Date<-substr(rownames(intdf),7,9)

cols<-ncol(guldf)-2

#calc sums
mid_sum<-aggregate(middf[,1:cols],  list(middf$Date), mean)
gul_sum<-aggregate(guldf[,1:cols],  list(guldf$Date), mean)
atl_sum<-aggregate(atldf[,1:cols],  list(atldf$Date), mean)
int_sum<-aggregate(intdf[,1:cols],  list(intdf$Date), mean)

rownames(mid_sum)<-paste("Mid",mid_sum[,1],sep="_")
rownames(gul_sum)<-paste("Gul",gul_sum[,1],sep="_")
rownames(atl_sum)<-paste("Atl",atl_sum[,1],sep="_")
rownames(int_sum)<-paste("Int",int_sum[,1],sep="_")

#combine
tot_sum<-rbind(mid_sum, gul_sum, atl_sum, int_sum)
tot_sum$Group.1<-NULL

#calc dists
tot_sam_rel<-decostand((tot_sum), method="hellinger")
tot_sam_bray<-vegdist(tot_sam_rel, method = "bray")
hist(tot_sam_bray)

pop_sites = c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4)
date_surveys=c(rep(1:5,4))

trajectoryPCoA(tot_sam_bray, pop_sites, date_surveys, traj.colors = c(4,2,1,3), lwd = 1)
legend("topleft", col=c(4,2,1,3), 
       legend=c("Midwest", "Gulf", "Atlantic", 
                "Intermediate"), bty="n", lty=1, lwd = 2)

#Trajectory stats####
#Are populations different in trajectory?

trajectoryDistances(tot_sam_bray, pop_sites, date_surveys, distance.type = "DSPD")
#          1         2         3
#2 0.2979923                    
#3 0.2926725 0.2246748          
#4 0.3302351 0.2330609 0.2568356

#trajectoryConvergence(tot_sam_bray, pop_sites, date_surveys)

#can we tell anything about the otus contributing to trajectory?
#remove otus sequentially, test how much they change the traj. 
#looping
drop_otu_tcomp<-function(tot_sums, otunum){
  tot_sum_1<-tot_sums
  tot_sum_1[,otunum]<-0
  
  pop_sites_1 = c(rep(1,5),rep(2,5),rep(3,5),rep(4,5),rep(5,5),rep(6,5),rep(7,5),rep(8,5))
  date_surveys_1=c(rep(1:5,8))
  
  comb_sum<-rbind(tot_sums, tot_sum_1)
  
  #calc dists
   tot_sam_rel_1  <-decostand((comb_sum), method="hellinger")
  tot_sam_bray_1  <-vegdist(tot_sam_rel_1, method = "bray")
  
  tdists<-trajectoryDistances(tot_sam_bray_1, pop_sites_1, date_surveys_1, distance.type = "DSPD")
  #out 1-4,2-5,3-6
  return(tdists[c(3,8,12)])
}

#run on top 2000 otus, nothing below this is prevalent enough to change trajectories
first_2000<-vector('list', 2000)
for(i in 1:2000){
  first_2000[[i]]<-drop_otu_tcomp(tot_sum, i)
}
f2000_means<-data.frame(otu_num=colnames(tot_sum)[1:2000], 
                       tdist_mean=sapply(first_2000, mean), 
                       tdist_max=sapply(first_2000, max))

top20<-sort(f2000_means$tdist_mean, decreasing = T)[20]

taxa_df<-tax_table(OTU_funkb)

taxa_df[rownames(taxa_df) %in% f2000_means$otu_num[f2000_means$tdist_mean>top20],]
trajectory_OTUs<-taxa_df[rownames(taxa_df) %in% f2000_means$otu_num[f2000_means$tdist_mean>top20],]
save(trajectory_OTUs, file="rdata/trajectory_OTUs_feb1.rda")

#Anova for traj####

####
# Average trajectory lengths by habitat
#### from Miquel de Caceres Ainsa
#load("Rdata/bci_ba_20_delta_bray.Rdata")
#need to get individual trajectories for each plant. subset to plant that have points repeated across all time points. 
load("~/Desktop/phyllos_all/rdata/dists_heatmap_jan27.rda")
repped<-rownames(dists$Time1_158)
meta$Code<-rownames(meta)
meta_repped<-meta[meta$PLANT_ID %in% repped,]
OTU_repped<-prune_samples(meta$PLANT_ID %in% repped,OTU_funkb)

#remove replicates (infected versus uninfected)
meta_repped$plotdate<-substr(meta_repped$Code,1,9)
meta_repped2<-meta_repped[!duplicated(meta_repped$plotdate), ]
OTU_repped2<-prune_samples(meta$Code %in% meta_repped2$Code,OTU_funkb)

repped_trans  <-decostand(t(otu_table(OTU_repped2)), method="hellinger")
repped_bray  <-vegdist(repped_trans, method = "bray", na.rm = T)
meta_repped2$Subpopulations<-as.factor(meta_repped2$Subpopulations)

repped_sites<-unique(meta_repped2$PLOT_GL)
repped_cols<-data.frame(sites=repped_sites,
                        pops=1)
repped_cols$pops<-meta_repped2$Subpopulations[match(repped_cols$sites,meta_repped2$PLOT_GL)]
repped_cols$pops<-as.factor(repped_cols$pops)
#Use modified trajectoryPCoA scripts from vegclust to allow for transparency and overplotting mean trajectories

trajectoryPCoA_mod <- function (d, sites, surveys = NULL, selection = NULL, traj.colors = NULL, 
                                axes = c(1, 2), asp,...) {
  siteIDs = unique(sites)
  nsite = length(siteIDs)
  if (is.null(selection)) 
    selection = 1:nsite
  else {
    if (is.character(selection)) 
      selection = (siteIDs %in% selection)
  }
  selIDs = siteIDs[selection]
  D2 = as.dist(as.matrix(d)[sites %in% selIDs, sites %in% selIDs])
  cmd_D2 <- cmdscale(D2, eig = TRUE, add = TRUE, k = nrow(as.matrix(D2)) - 
                       1)
  x <- (cmd_D2$points[, axes[1]])
  y <- cmd_D2$points[, axes[2]]
  plot(-x, y, type = "n", asp = asp, xlab = paste0("PCoA ", axes[1], 
                                                   " (", round(100 * cmd_D2$eig[axes[1]]/sum(cmd_D2$eig)), 
                                                   "%)"), ylab = paste0("PCoA ", axes[2], " (", round(100 * 
                                                                                                        cmd_D2$eig[axes[2]]/sum(cmd_D2$eig)), "%)"))
  sitesred = sites[sites %in% selIDs]
  if (!is.null(surveys)) 
    surveysred = surveys[sites %in% selIDs]
  else surveysred = NULL
  for (i in 1:length(selIDs)) {
    ind_surv = which(sitesred == selIDs[i])
    if (!is.null(surveysred)) 
      ind_surv = ind_surv[order(surveysred[sitesred == 
                                             selIDs[i]])]
    for (t in 1:(length(ind_surv) - 1)) {
      niini = ind_surv[t]
      nifin = ind_surv[t + 1]
      if (!is.null(traj.colors)) 
        arrows(x0=-(x[niini]), y[niini], x1=-(x[nifin]), y[nifin], 
               col = traj.colors[i], ...)
      else arrows(x[niini], y[niini], x[nifin], y[nifin], 
                  ...)
    }
  }
  invisible(cmd_D2)
}
trajectoryPCoA_mod2<-function (d, sites, surveys = NULL, selection = NULL, traj.colors = NULL, 
                               axes = c(1, 2), ...) {
  siteIDs = unique(sites)
  nsite = length(siteIDs)
  if (is.null(selection)) 
    selection = 1:nsite
  else {
    if (is.character(selection)) 
      selection = (siteIDs %in% selection)
  }
  selIDs = siteIDs[selection]
  D2 = as.dist(as.matrix(d)[sites %in% selIDs, sites %in% selIDs])
  cmd_D2 <- cmdscale(D2, eig = TRUE, add = TRUE, k = nrow(as.matrix(D2)) - 
                       1)
  x <- (cmd_D2$points[, axes[1]])
  y <- cmd_D2$points[, axes[2]]
  #plot(-x, y, type = "n", asp = 1, xlab = paste0("PCoA ", axes[1], 
  #                                               " (", round(100 * cmd_D2$eig[axes[1]]/sum(cmd_D2$eig)), 
  #                                               "%)"), ylab = paste0("PCoA ", axes[2], " (", round(100 * 
  #                                                                                                    cmd_D2$eig[axes[2]]/sum(cmd_D2$eig)), "%)"))
  sitesred = sites[sites %in% selIDs]
  if (!is.null(surveys)) 
    surveysred = surveys[sites %in% selIDs]
  else surveysred = NULL
  for (i in 1:length(selIDs)) {
    ind_surv = which(sitesred == selIDs[i])
    if (!is.null(surveysred)) 
      ind_surv = ind_surv[order(surveysred[sitesred == 
                                             selIDs[i]])]
    for (t in 1:(length(ind_surv) - 1)) {
      niini = ind_surv[t]
      nifin = ind_surv[t + 1]
      if (!is.null(traj.colors)) 
        arrows(x0=(x[niini]), y[niini], x1=(x[nifin]), y[nifin], 
               col = traj.colors[i], main="b")
      else arrows(x[niini], y[niini], x[nifin], y[nifin], 
                  ...)
    }
  }
  invisible(cmd_D2)
}

library(colorspace)
c1<-adjust_transparency(1, alpha=.2)
c2<-adjust_transparency(2, alpha=.2)
c3<-adjust_transparency(3, alpha=.2)
c4<-adjust_transparency(4, alpha=.2)

repped_cols$cols<-c(c1,c2,c3,c4)[match(repped_cols$pops, c("Atlantic", "Gulf", "Intermediate", "Midwest"))]

png(filename = "traj_arrows_jan28.png", width = 10,height = 6,units = "in",res = 400)
trajectoryPCoA_mod(repped_bray, meta_repped2$PLOT_GL, meta_repped2$Date, lwd = 1,length=.1, 
                   traj.colors = repped_cols$cols)

trajectoryPCoA_mod2(tot_sam_bray, pop_sites, date_surveys, traj.colors = c(4,2,1, 3), lwd = 2, length=.2)
legend("bottomright", col=1:4, 
       legend=c( "Atlantic", "Gulf","Intermediate", "Midwest"), bty="n", lty=1, lwd = 2)
graphics.off()

pcoa2<-cmdscale(repped_bray,k = 4)
plot(-pcoa2[,1],pcoa2[,2], col=(meta_repped2$Subpopulations))

#get lengths and directionality
quadrat_tl = trajectoryLengths(repped_bray, sites=meta_repped2$PLOT_GL, surveys = meta_repped2$Date)
quadrat_td = trajectoryDirectionality(repped_bray, sites = meta_repped2$PLOT_GL, surveys=meta_repped2$Date)

#get mean trajectory lengths by population and evaluate anova by factor
#check order is the same
rownames(quadrat_tl)==repped_cols$sites
tapply(quadrat_tl$Trajectory, repped_cols$pops, mean, na.rm=T)
mod= aov(quadrat_tl$Trajectory~repped_cols$pops)
summary(mod)
TukeyHSD(mod)
plot(TukeyHSD(mod))

#get mean traj direction and evaluate aov
tapply(quadrat_td, repped_cols$pops, mean, na.rm=T)
mod2= aov(quadrat_td~repped_cols$pops)
summary(mod2)
TukeyHSD(mod2)
plot(TukeyHSD(mod2))

#MTV-LMM ####
#https://github.com/cozygene/MTV-LMM

#set up files
#count table - A matrix of temporal samples by taxa, across multiple hosts. The 
#first row contains the sample ids. The 
#first column includes taxa ids. Then every consecutive column contains 
#read counts for each sample. Note that this order must be respected (see example below).
#from A.
load("~/Desktop/phyllos_all/rdata/OTU_funkb_jan31.rda")

meta1<-sample_data(OTU_funkb)
meta2<-as.data.frame(do.call(cbind, meta1@.Data))
colnames(meta2)<-meta1@names

#need to trim down to ones that are in all time points and non-duplicated
#split by time
plant_idlist<-split(meta2, as.character(meta2$Date))

shared<-intersect(plant_idlist$`158`$PLANT_ID,
                  intersect(plant_idlist$`212`$PLANT_ID,
                  intersect(plant_idlist$`233`$PLANT_ID, 
                  intersect(plant_idlist$`260`$PLANT_ID, 
                  plant_idlist$`286`$PLANT_ID))))

OTU_shared<-prune_samples(meta2$PLANT_ID %in% shared,OTU_funkb)

#trim duplicated

meta3<- meta2[meta2$PLANT_ID %in% shared,]
meta3$PD<-paste(meta3$PLOT_GL,meta3$Date, sep="_")

OTU_unique<-prune_samples(!duplicated(meta3$PD), OTU_shared)

count_table1<-otu_table(OTU_unique)
count_table2<-as.data.frame(count_table1)
count_table2[1:5,1:5]

#metadata - The 
#first row contains the headers ('sample_id', 'ind_id', 'Id_num', 'ind_time', 'Sampling_day'). The 
#first column contains the sample ids. The 
#second column contains the subject ids; the 
#third column is the index of each subject (between 1 - number of subjects in the data). The 
#fourth column is the time index where the scale is the experiment's sampling rate (e.g., days, weeks, months). The 
#fifth column is the sampling day (similar to the fourth column, but on the scale of days). Note that these names must be respected (see examples below).

metax1<-sample_data(OTU_unique)
metax2<-as.data.frame(do.call(cbind, metax1@.Data))
metax3<-metax2
colnames(metax3)<-metax1@names
metax3$sample_id<-metax1@row.names
metax3$Id_num<-as.factor(metax3$PLOT_GL)
sort(unique(as.numeric(metax3$Id_num)))
metax3$Id_num<-as.numeric(metax3$Id_num)
metax3$Sampling_day<-as.numeric(metax3$Date)
#ind time seems to be just scaled by month
metax3$ind_time<-metax3$Sampling_day/30
table(metax3$ind_time)

#rearrange cols
metax4<-metax3[,c("sample_id", "PLOT_GL", "Id_num", "ind_time", "Sampling_day")]
colnames(metax4)[2]<-"ind_id"
rownames(metax4)[duplicated(rownames(metax4))]
table(metax4$ind_id)
write.csv(metax4, file = "rdata/meta4_feb1.csv",row.names = F)

#in bash ./run_TE.sh ~/Downloads/gcta_1.93.2beta_mac/gcta64 ~/Desktop/phyllos_all/raw_outs/count_table_test.csv ~/Desktop/phyllos_all/raw_outs/meta4.csv 

results1<-read.csv("~/MTV-LMM/MTV_LMM/Results/TE_Results.csv")
taxa<-rownames(count_table3)[results1$taxa_index]
results1$otus<-taxa

res_sort<-results1[order(results1$Time_explainability, decreasing = T),]
res_sort$otus[1:20]

res_sort2<-results1[order(results1$logL, decreasing = T),]
res_sort2$otus[1:20]


plot(1:nrow(res_sort), res_sort$Time_explainability)

bonnf<-0.05/(nrow(results1))

res_top<-res_sort[which(res_sort$p_value_adjusted<bonnf),]

write_csv(res_top, "rdata/MTV-LMM_top_feb1.csv")


#Venn####
#devtools::install_github("yanlinlin82/ggvenn")
library(ggvenn)

MTV_otus<-MTV_LMM_top_feb1$otus
traj_otus<-rownames(trajectory_OTUs@.Data)

core_taxa <- read_csv("~/Downloads/core_taxa.csv")
core_otus<-core_taxa$x

otusums<-rowSums(otu_table(OTU_funkb))
topsums<-names(sort(otusums, decreasing = T)[1:100])

over1<-intersect(MTV_otus, traj_otus)
overall<-intersect(MTV_otus, intersect(core_otus,traj_otus))

venn1<-list(MTV_LMM=MTV_otus,
            Trajectory=traj_otus,
            #Abundant=topsums,
            Core=core_otus)

venp<-ggvenn(venn1,
       fill_color = c("skyblue", "orange", "tomato", "pink"),
       stroke_size = .5,text_size = 7,show_percentage = F)

venp + plot_annotation(tag_levels = 'a')
ggsave("../figs/venn_plot_feb28.png", height = 4, width=4)
ggsave("../figs/venn_plot_feb28.pdf", height = 4, width=4)
ggsave("../figs/venn_plot_feb28.eps", height = 4, width=4)

#core taxa####
core_tax<-prune_taxa(rownames(OTU_funkb@tax_table) %in% core_otus, OTU_funkb)
View(core_tax@tax_table@.Data)

core_taxa_ids<-as.data.frame(core_tax@tax_table@.Data)

write_csv(core_taxa_ids, "rdata/core_taxa_ids.csv")

otutab_core<-otu_table(core_tax)

write.csv(otutab_core,file =  "otutab_core.csv")


#devtools::install_github("ropenscilabs/datastorr")
#devtools::install_github("traitecoevo/fungaltraits")
library(fungaltraits)

traits<-fungal_traits()

#sub to genera in core

traits_core<-traits[traits$Genus %in% core_taxa_ids$Genus,]

traits_core_trim<-traits_core %>% select(Genus,species, growth_form_fg, guild_fg, notes_fg, source_funguild_fg)

#try MTV taxa
MTV_tax<-prune_taxa(rownames(OTU_funkb@tax_table) %in% MTV_otus, OTU_funkb)
MTV_taxa_ids<-as.data.frame(MTV_tax@tax_table@.Data)
traits_MTV<-traits[traits$Genus %in% MTV_taxa_ids$Genus,]
traits_MTV_trim<-traits_MTV %>% select(Genus,species, growth_form_fg, guild_fg, notes_fg, source_funguild_fg)

#overlap
overall_tax<-prune_taxa(rownames(OTU_funkb@tax_table) %in% overall, OTU_funkb)
overall_taxa_ids<-as.data.frame(overall_tax@tax_table@.Data)
traits_overall<-traits[traits$Genus %in% overall_taxa_ids$Genus,]
traits_overall_trim<-traits_overall %>% select(Genus,species, growth_form_fg, guild_fg, notes_fg, source_funguild_fg)

overall_write<-overall_taxa_ids %>% select(Genus, BestMatch, HL_hit_percent_id)%>%arrange(Genus)%>%mutate(HL_hit_percent_id=round(as.numeric(HL_hit_percent_id), 1))
colnames(overall_write)<-c("Genus","Best_Match", "BLAST_percent_ID")
overall_write$Guild_Estimate<-c("Pathogen", "Mycoparasite","Mycoparasite", "Pathogen", "Pathogen", "Pathogen", "Mycoparasite", "Yeast", "Yeast", "Yeast", "Yeast", "Pathogen", "unknown", "unknown")

#library(kableExtra)

overall_write %>%
  kbl(caption = "b") %>%
  kable_classic(full_width = F, html_font = "Cambria")%>%
  #row_spec((1:(nrow(adonis_tests2[[4]])-2)),bold=T,hline_after = F)%>%
  save_kable("../figs/OTU_overlap_table_mar23.pdf")

overall_write %>%
  kbl(caption = "b") %>%
  kable_classic(full_width = F, html_font = "Cambria")%>%
  #row_spec((1:(nrow(adonis_tests2[[4]])-2)),bold=T,hline_after = F)%>%
  save_kable("../figs/OTU_overlap_table_mar23.html")

#add in indicator species
indicators<-c("OTU_6", "OTU_52", "OTU_4", "OTU_1723", "OTU_15", "OTU_1207")

ind_tax<-tax_all[tax_all$OTU_ID %in% indicators,]

colnames(overall_write)
core_taxa_ids <- read.csv("~/Downloads/core_taxa_ids - core_taxa_ids (1).csv")

ind_tax$Guild_Estimate<-core_taxa_ids$functional_grp[match(ind_tax$OTU_ID, core_taxa_ids$OTU_ID)]

ind_tax2 <- ind_tax %>%
  select(Genus, BestMatch, HL_hit_percent_id, Guild_Estimate) %>%
  mutate(HL_hit_percent_id=round(as.numeric(HL_hit_percent_id), 1))%>%
  arrange(BestMatch)

(ind_tax2)[rownames(ind_tax2) %in% c("OTU_4", "OTU_6"),]
rownames(ind_tax2)[rownames(ind_tax2) %in% c("OTU_4", "OTU_6")]<-c("OTU_6.", "OTU_4.")
#check order
(ind_tax2)[rownames(ind_tax2) %in% c("OTU_4.", "OTU_6."),]


colnames(ind_tax2)<-colnames(overall_write)
#add spacing col.
spacer<-ind_tax2[1,]
spacer[1,]<-" "
rownames(spacer)<- "Infection Indicator Taxa"

write2<-rbind(overall_write, spacer,ind_tax2)

write2 %>%
  kbl(caption = "b") %>%
  kable_classic(full_width = F, html_font = "Cambria")%>%
  row_spec(15,bold=T,hline_after = F)%>%
  save_kable("figs/OTU_overlap_table_mar23.pdf")

write2$OTU_ID<-rownames(write2)

write_csv(write2, "write2_mar23.csv" )

#Rare taxa####
claviceps<-rownames(taxtab)[which(taxtab$Genus=="Claviceps")]
clavtab<-otu_table_all_lib[otu_table_all_lib[,1] %in% claviceps,]
clavsum<-colSums(clavtab)
clavsum[clavsum>0]

meta<-rownames(taxtab)[which(taxtab$Genus=="Metarhizium")]
metatab<-otu_table_all_lib[otu_table_all_lib[,1] %in% meta,]
metasum<-colSums(metatab[,-1])
metasum[metasum>0]
length(metasum[metasum>0])

coll<-rownames(taxtab)[which(taxtab$Genus=="Colletotrichum")]
colltab<-otu_table_all_lib[otu_table_all_lib[,1] %in% coll,]
collsum<-colSums(colltab[,-1])
collsum[collsum>0]
length( collsum[collsum>0])

save.image("~/Desktop/phyllos_all/rdata/F.phyllos_taxa_feb11.RData")
