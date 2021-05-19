 ##Phyllosphere clean species summaries
#Acer VanWallendael, cleaned to trim & non-fungi

library(readr)
library(tidyverse)
library(readxl)

setwd("~/Desktop/phyllos_all/rdata")

#Clean OTU tables####
otu_table_all_lib <- as.data.frame(read_delim("~/Desktop/phyllos_all/raw_outs/otu_table_all_lib_ITS_UPARSE_R1.txt",  "\t", escape_double = FALSE, trim_ws = TRUE))

#decontam####
library(decontam)
library(phyloseq)
library(ggplot2)

otu_table<-otu_table_all_lib

  #clean code
  otu_t<-t(otu_table)[-1,]
  colnames(otu_t)<-otu_table[,1]
  
  rownames(otu_table)<-otu_table[,1]
  otu_table<-otu_table[,-1]
  
  #need to fix some labels
  rownames(otu_t)[substr(rownames(otu_t),1,2)=="44"]<-"M4418x260I"
  rownames(otu_t)[substr(rownames(otu_t),7,8)=="77"]<-"K2104x077U"
  rownames(otu_t)[substr(rownames(otu_t),1,5)=="M5811"]<-"M5817x286U"
  rownames(otu_t)[substr(rownames(otu_t),1,10)=="M1918x077U"]<-"K1918x077U"
  rownames(otu_t)[substr(rownames(otu_t),1,5)=="M2617"]<-"M2607x233I" 
  rownames(otu_t)[substr(rownames(otu_t),1,5)=="M5419"]<-"M5418x260U"
  colnames(otu_table)<-rownames(otu_t)
  
  codes<-rownames(otu_t)
  sampledata<-data.frame(
      PLOT_GL   =substr(codes,1,5),
      Date      =as.numeric(substr(codes,7, 9)),
      site      =substr(codes,1,1),
      inf       =substr(codes,10,10),
      ident     =substr(codes,1,2)
  )
  rownames(sampledata)<-codes
  
  #read in metadata
  DOE <- read_excel("~/Downloads/DOE_GWAS_Master Plant List (2).xlsx")
  PVDIV <- read_excel("~/Downloads/PVDIV_Master Metadata File_5-8-2019.xlsx")
  
  sampledata$PLANT_ID<-DOE$PLANT_ID[match(sampledata$PLOT_GL, DOE$PLOT_GL)]
  sampledata$Eco<-PVDIV$SUBPOP_SNP[match(sampledata$PLANT_ID, PVDIV$PLANT_ID)]
  unique(sampledata$Eco)
  subpops<-unique(sampledata$Eco)
  bigpops<-c("Gulf", "Midwest", NA, "Atlantic", "Intermediate",  "Gulf", "Intermediate", "Intermediate", "Intermediate")
  sampledata$Subpopulations<-bigpops[match(sampledata$Eco, subpops)]
  
  ##set up phyloseq object
  
  OTU = otu_table(otu_table, taxa_are_rows = T)
  
  #need sample vars. rows are sample_names
  SAMPLES = sample_data(sampledata)
  
  
  #need phylo assignments for each OTU. 
  load("~/Downloads/constax_taxonomy_07.RData")
  #use new taxonomy. check if IDs match up
  constax_taxonomy_07$OTU_ID[1:50] == rownames(taxmat)[1:50]
  constax_taxonomy_07$Genus[1:20] == taxmat[,6][1:20]
  
  TAX<-tax_table(as.matrix(constax_taxonomy_07))
  
  phylo_obj<-phyloseq(OTU, SAMPLES, TAX)
  
  ps <- phylo_obj
  #BLANK118T is not a true blank 
  ps <- prune_samples(ps@sam_data@row.names!="BLANK118T", ps)
  #"M6610x159U" "M3803x159U" "M2614x159U" are mislabeled
  ps <- prune_samples(ps@sam_data@row.names!="M6610x159U", ps)
  ps <- prune_samples(ps@sam_data@row.names!="M3803x159U", ps)
  ps <- prune_samples(ps@sam_data@row.names!="M2614x159U", ps)
  
  df1 <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
  df1$LibrarySize <- sample_sums(ps)
  #take out replicated library preps with lower seq depth. 
  df1$PLOT_DATE<-substr(row.names(df1),1,10)
  reptab<-table(df1$PLOT_DATE)
  reps<-names(reptab[reptab>1])
  repdf1<-as.data.frame(do.call(cbind,df1[df1$PLOT_DATE %in% reps,]@.Data))
  colnames(repdf1)<-df1[df1$PLOT_DATE %in% reps,]@names
  rownames(repdf1)<-df1[df1$PLOT_DATE %in% reps,]@row.names
  repdf1$LibrarySize<-as.numeric(repdf1$LibrarySize)
  replist<-split(repdf1, repdf1$PLOT_DATE)
  reprm<-vector("list", length(replist))
  sorted<-vector("list", length(replist))
  for(i in 1:length(replist)){
    sorted<-replist[[i]][order(replist[[i]]$LibrarySize, decreasing = T),]
    reprm[[i]]<-rownames(sorted[-1,])
  }
  reprms<-do.call(c, reprm)
  ps.reprm<-prune_samples(!(rownames(df1) %in% reprms),ps)
  
  
#DECONTAMINATION####
  
  df <- as.data.frame(sample_data(ps.reprm)) # Put sample_data into a ggplot-friendly data.frame
  df$LibrarySize <- sample_sums(ps.reprm)
  df$orig_order<-c(1:nrow(df))
  #REORDERING STEP
  df <- df[order(df$LibrarySize),]
  df$Index <- seq(nrow(df))
  df$control<-df$ident
  df$control[df$control %in% c("EM","BL")]<-"Control"
  df$control[df$control %in% c("NE","Ne")]<-"Negative"
  df$control[df$control == "Mo"]<-"Mock"
  df$control[df$control %in% c("M6",      "M3",      "M5" ,     "C6" ,     "B2"  ,    "K1" ,     "M1" ,    
                               "M2",      "P5",      "C4" ,     "K2" ,     "M4"  ,    "C7" ,     "P6" ,     "B1" ,    
                               "C5",      "P2",      "P7" ,     "B3" ,     "P4"  ,    "C2" ,     "P3" ,     "P1"     
                                )]<-"Samp"
  df <- df[order(df$orig_order),]
  ggplot(data=df, aes(x=Index, y=LibrarySize, color=control)) + 
    geom_jitter(width=2, shape=1)+
    scale_color_manual(values = c("red", "black", "purple","yellow"))+
    theme_classic()
  #ggsave("../figs/indexed_contams_jan27.png", height = 5, width=5)
  #ggsave("../figs/indexed_contams_jan27.pdf", height = 5, width=5)
  #ggsave("../figs/indexed_contams_jan27.eps", height = 5, width=5)
  
  #logical of controls
  sample_data(ps.reprm)$is.neg <- df$control %in% c("Control", "Negative")
  #contamdf.prev <- isContaminant(ps.reprm, method="prevalence", neg="is.neg")
  #table(contamdf.prev$contaminant)
 #head(which(contamdf.prev$contaminant))
  
  #raise the threshold
  contamdf.prev05 <- isContaminant(ps.reprm, method="prevalence", neg="is.neg", threshold=0.5)
  table(contamdf.prev05$contaminant)
  
  #82 contaminants, 7881 otus
  
  # Make phyloseq object of presence-absence in negative controls and true samples
  ps.pa <- transform_sample_counts(ps.reprm, function(abund) 1*(abund>0))
  ps.pa.neg <- prune_samples(df$control %in% c("Control","Negative"), ps.pa)
  ps.pa.pos <- prune_samples(df$control %in% c("Samp","Mock"), ps.pa)
  # Make data.frame of prevalence in positive and negative samples
  df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                      contaminant=contamdf.prev05$contaminant)
  ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + 
    geom_point(size=1,  alpha=.5 ) + 
    xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")+
    scale_color_manual(values = c("skyblue", "firebrick"))+
    theme_classic()
  #ggsave("../figs/prev_contams_jan27.png", height = 4, width=5)
  #ggsave("../figs/prev_contams_jan27.pdf", height = 4, width=5)
  #ggsave("../figs/prev_contams_jan27.eps", height = 4, width=5)

  
#remove contaminants with the higher threshold
decont_all_lib<-prune_taxa(contamdf.prev05$contaminant==FALSE, ps.reprm)

taxa_df<-as.data.frame(tax_table(decont_all_lib))
#how many levels otus of each were IDed
apply(taxa_df, 2, FUN = function(x){length(x[!is.na(x)])})

#OTU_ID             Kingdom              Phylum               Class 
#7881                7414                6870                6457 
#Order              Family               Genus             Species 
#6171                5106                3658                2131 

#Non-normalized phyloseqs####
ps.samp <- prune_samples(!(decont_all_lib@sam_data@.Data[[5]] %in% c("Ne","NE", "BL", "EM","Mo")), decont_all_lib)

#fungi_only
fungi_nonnorm<-prune_taxa(ps.samp@tax_table@.Data[,2]=="Fungi", ps.samp)
KBS_fungi_nonnorm<-prune_samples(fungi_nonnorm@sam_data@.Data[[3]]=="M",fungi_nonnorm)

#Normalize####
#from https://github.com/Gian77/Scientific-Papers-R-Code/blob/410fac45037da835d23dc906c703a2763ef761b0/Benucci_etal_2020_Patient_Propagules/S2_File_R_scripts_revised.R

CSSNorm <-function(dataframe){
  require(metagenomeSeq)
  dataframe %>% 
    phyloseq_to_metagenomeSeq() -> physeq_CSS
  p_biom <-cumNormStatFast(physeq_CSS)
  biom_quant <-cumNorm(physeq_CSS, p=p_biom)
  physeq_CSS <- MRcounts(biom_quant, norm=T)
  physeq_mSeq <- dataframe
  otu_table(physeq_mSeq) <- otu_table(physeq_CSS, taxa_are_rows=TRUE)
  return(physeq_mSeq)
}

#warning: sample with few features

#remove controls
ps.samp <- prune_samples(!(decont_all_lib@sam_data@.Data[[5]] %in% c("Ne","NE", "BL", "EM","Mo")), decont_all_lib)
#normalize
OTU_scale<-CSSNorm(ps.samp)
OTU_shared<-CSSNorm(fungi_shared)

#check. the vals should now be decimal
#hist(log(apply(ps.samp@otu_table@.Data,2,mean)), main="Pre-Norm")
#hist(log(apply(OTU_scale@otu_table@.Data,2,mean)), main="Post-norm")

#ps.samp@otu_table@.Data[1:5,1:5]
#OTU_scale@otu_table@.Data[1:5,1:5]

#Check samples####
#totcov<-apply(OTU_scale@otu_table@.Data,2,sum)
#sort(totcov, decreasing = T)[1:50]

#make an object with only fungi
OTU_fungi<-prune_taxa(OTU_trim1@tax_table@.Data[,2]=="Fungi",OTU_trim1)

#make an object with only focal site
OTU_funkb<-prune_samples(OTU_fungi@sam_data@.Data[[3]]=="M", OTU_fungi)

##get shared across sites #####
#here are the IDs of the samples shared across sites
allsites<- c("C7020", "C6911", "M1202", "K2819", "P7020", "C7011", "P2614", "M1801", "M4211", "M5405", "C4806", "P7413", "K1204", "K1918", "P2806", "K1414")
#pull out all +kbs eqs
meta_sampall <- as.data.frame(do.call(cbind,fungi_nonnorm@sam_data@.Data))
colnames(meta_sampall)<-c("PLOT_GL", "Date","Site", "inf","ID", "PLANT_ID", "Eco", "Subpopulation", "Contam")
row.names(meta_sampall)<-(fungi_nonnorm@sam_data@row.names)

fungi_shared1<-prune_samples(meta_sampall$PLOT_GL %in% allsites, fungi_nonnorm)

#take out Brookings (not enough samples)
fungi_shared<-prune_samples(fungi_shared1@sam_data@.Data[[3]] != "B", fungi_shared1)
#Outputs####
save(fungi_nonnorm,      file = "fungi_nonnorm_jan31.rda")
save(KBS_fungi_nonnorm,file="KBS_fungi_nonnorm_feb18.rda")
save(decont_all_lib,      file="decont_all_lib_jan31.rda")  
save(OTU_scale,                file="OTU_scale_jan31.rda")
save(OTU_fungi,                file="OTU_fungi_jan31.rda")
save(OTU_funkb,                file="OTU_funkb_jan31.rda")
save(OTU_shared,              file="OTU_shared_jan31.rda")
#save.image(         "A.phyllos_clean_summaries.RData")


#sample table for SRA#####
namesall<-(rownames(t(otu_table_all_lib))
            )
names<-data.frame(samp=namesall, title=substr(namesall,1,4))

names2<-names[-1,]
names2$bioproj<-NA
names2$Organism<-"Panicum virgatum"
names2$Collection_date<-substr(names2$samp,7,9)
names2$Collection_date[(nchar(names2$Collection_date))<3]<-NA
names2$context<-"Leaf Tissue"
names2$context[is.na(names2$Collection_date)]<-"Control"
names2$context_local<-names2$context
names2$context<-"Field"
names2$env_medium<-names2$context
names2$env_medium<-names2$context_local
names2$location<-substr(names2$samp,1,1)
unique(names2$location)
names2$location<-c("Columbia MO", "KBS MI", "Control", "Kingsville TX","Austin TX", "Brookings SD", "KBS MI", "Control")[match(names2$location, unique(names2$location))]
names2$host<-"Panicum virgatum"
names2$longlat<-names2$location
names2$longlat[names2$location=="KBS MI"]<-"42.419, -85.371"
names2$longlat[names2$location=="Columbia MO"]<-"38.896, -92.217"
names2$longlat[names2$location=="Austin TX"]<-"30.383, -97.729"
names2$longlat[names2$location=="Kingsville TX"]<-"27.549, -97.881"
names2$longlat[names2$location=="Brookings SD"]<-"44.352, -96.803"
names2$longlat[names2$location=="Control"]<-NA
write_csv(names2, "names2_ncbi_sra_mar17.csv")
