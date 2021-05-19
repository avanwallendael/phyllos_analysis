#RNA Expression data for GWAS outlier region
#Acer VanWallendael

setwd("~/Desktop/GWAS")

meta1 <- read.csv("~/Downloads/Archive/md4geno3site4Acer.csv")
counts1 <- read.csv("~/Downloads/Archive/cnts4geno3site4Acer.csv")

library(DESeq2)
library(tidyverse)

#Outlier-adjacent genes
#Pavir.2NG521906.1, Pavir.2NG521912.1, Pavir.2NG521915.1

counts1[1:5,1:5]

#clean
foc<-counts1[c("Pavir.2NG521906", "Pavir.2NG521912", "Pavir.2NG521915"),]
foc$id<-rownames(foc)

foc_long<-pivot_longer(foc, -id)
foc_long$library<-foc_long$name
foc_long2<-left_join(foc_long, meta1, by="library")

#test plot
ggplot(foc_long2, aes(x=site, y=value, col=geno))+
  geom_jitter()+
  facet_grid(.~id)+
  theme_classic()

foc_long2$eco<-as.character(foc_long2$geno)
foc_long2$eco[foc_long2$eco %in% c("AP13", "WBC")]<-"Lowland"
foc_long2$eco[foc_long2$eco %in% c("DAC", "VS16")]<-"Upland"

table(foc_long2$eco)

foc_long2$geno<-factor(foc_long2$geno, levels=c("AP13",'WBC',"DAC", "VS16"))
foc_long2$Site<-foc_long2$site
foc_long2$Site[foc_long2$Site=="PKLE"]<-"Pickle"
foc_long2$Site[foc_long2$Site=="CLMB"]<-"Columbia"
foc_long2$Site[foc_long2$Site=="KBSM"]<-"KBS"

foc_long2$Site<-factor(foc_long2$Site, levels = c("KBS","Columbia", "Pickle"))


ggplot(foc_long2, aes(x=geno, y=value, fill=eco))+
  geom_boxplot(fill=NA)+
  geom_jitter(width=.2,  aes(col=eco))+
  scale_fill_manual(values=c("firebrick", "skyblue"), name="Ecotype")+
  scale_color_manual(values=c("firebrick", "skyblue"), name="Ecotype")+
  labs(x="Genotype", y="Transcript Counts")+
  facet_grid(id~Site, scales = "free")+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle = 90))

ggsave("figs/rna_boxes2_feb15.png", height = 4,width=6)
ggsave("figs/rna_boxes2_feb15.pdf", height = 4,width=6)
ggsave("figs/rna_boxes2_feb15.eps", height = 4,width=6)


#Deseq without datatable####

#trim to focal
foc1<-counts1[c("Pavir.2NG521906", "Pavir.2NG521912", "Pavir.2NG521915"),]
meta1$eco<-meta1$geno
meta1$eco[meta1$eco %in% c("AP13", "WBC")]<-"Lowland"
meta1$eco[meta1$eco %in% c("DAC", "VS16")]<-"Upland"

#drop low counts bc they can impact things
foc2<-foc1
foc2[1,][foc2[1,]>0&foc2[1,]<6]<-5
foc2[2,][foc2[2,]>0&foc2[2,]<6]<-5
foc2[3,][foc2[3,]>0&foc2[3,]<6]<-5

#obj construction
ds_test<-DESeqDataSetFromMatrix(foc2, meta1, design = ~eco)

wd.grp <- DESeq(
  #object DESeqDataSet DESeqDataSetFromHTSeqCount
  object=ds_test,
  #"Wald" or "LRT
  test="Wald",
  #design. don't need to specify if you add to deseq obj
 # full = ~ eco,
  #null design only for LRT
  #reduced= ~ 1,
  #print?
  quiet = F,
  parallel = F)

results(wd.grp)

meta1$ecosite<-paste0(meta1$eco,meta1$site)
#https://biohpc.cornell.edu/doc/RNA-Seq-2019-Lecture3.pdf

ds_test2<-DESeqDataSetFromMatrix(foc2, meta1, design = ~ eco + site +eco:site)

lrt_site<-DESeq(object = ds_test2, test="LRT", reduced = ~eco + site, quiet = F)
results(lrt_site)


res1<-results(wd.grp, format = "DataFrame")
res2<-do.call(cbind,res1@listData)
rownames(res2)<-res1@rownames

library(kableExtra)
#reformat
res2<-as.data.frame(res2)
res3<-res2 %>%
  mutate(across(1:4,round, 2))
  
colnames(res3)<-c("Base Mean", "log Fold Change", "SE", "W", "p-value", "Adj. p-value")
  
  res3 %>% kbl(caption = "Wald Test") %>%
  kable_classic(full_width = F, html_font = "Cambria")%>%
  #row_spec((1:(nrow(adonis_tests2[[4]])-2)),bold=T,hline_after = F)%>%
  save_kable("figs/wald_feb10.pdf")