##Phyllosphere gwas
#Acer VanWallendael, Nov 12 2020

setwd("~/Desktop/phyllos_all/rdata")

library(switchgrassGWAS)
library(bigsnpr)
library(tidyverse)

#load snps
snp <- snp_attach("~/Downloads/Pvirgatum_V5_GWAS_630g_33M_tensite_twoyear.rds")
CHRN<-snp$map$chromosome
chrkey<-data.frame(chr=unique((CHRN)), 
                   num=1:length(unique(CHRN)))
chrnum<-chrkey$num[match(CHRN, chrkey$chr)]
POS<-snp$map$physical.pos

#load distmat
load("~/Desktop/GWAS/rdata/gwas_mat_all.rda")
gwas_mat_all<-as.matrix(gwas_mat_all[,1:785])
#correct for population structure with Single-Vector decomposition of snps
#need to convert CHROMs to numeric
CHRN<-snp$map$chromosome
chrkey<-data.frame(chr=unique((CHRN)), 
                   num=1:length(unique(CHRN)))
chrnum<-chrkey$num[match(CHRN, chrkey$chr)]
POS<-snp$map$physical.pos

obj.svd<-snp_autoSVD(snp$genotypes, infos.chr = chrnum, infos.pos = POS)

#load phenos
load("~/Desktop/phyllos_all/rdata/points3_otus_jan27.rda")

#sub to one time pt with max genetic separation
phenos_gwas1<-points3[which(points3$Date==260),]
#rearrange
phenos_gwas2<-as.data.frame(phenos_gwas1[,c("PLANT_ID", "MDS2")])
#split to all PLANT_IDs
all_IDs<-data.frame(PLANT_ID=snp$fam$sample.ID)
phenos_gwas3<- left_join(all_IDs, phenos_gwas2, by="PLANT_ID")

#sub mat to plant ids
#check <- attr(gwas_mat_all, "dimnames")[[2]] %in% phenos_gwas2$PLANT_ID

phyllos_gwas<-pvdiv_gwas(phenos_gwas3, type = "linear",snp = snp,
                         #ncores = nb_cores(),
                         covar = obj.svd)

#predict the pvals
pvals<-stats::predict(phyllos_gwas)

#replace the NAs
pvals[is.na(pvals)]<-1

phyllos_gwas$pvals<-pvals
phyllos_gwas$log10p<-(-(pvals))
phyllos_gwas$p<- 10^pvals

phyllos_gwas$FDR_adj <- p.adjust(phyllos_gwas$p, method = "BH")

FDRthreshhi <- phyllos_gwas %>%
  as_tibble() %>%
  #for fdr between .1 and 1, what is the max log10p?
  #filter(between(.data$FDR_adj, 0.10001, 1)) %>%  # 10% FDR threshold
  filter(between(.data$FDR_adj, 0.05001, 1)) %>%  # 5% FDR threshold
  summarise(thresh = max(.data$log10p))
FDRthreshlo <- phyllos_gwas %>%
  as_tibble() %>%
  #for fdr between 0 and 0.0999, what is the min log10p?
  #filter(between(.data$FDR_adj, 0, 0.09999)) %>%  # 10% FDR threshold
  filter(between(.data$FDR_adj, 0, 0.04999)) %>%  # 5% FDR threshold
  summarise(thresh = min(.data$log10p))
if(FDRthreshhi$thresh[1] > 0 & FDRthreshlo$thresh[1] > 0){
  FDRthreshold = (FDRthreshhi$thresh[1] + FDRthreshlo$thresh[1])/2
} else if(FDRthreshhi$thresh[1] > 0){
  FDRthreshold = FDRthreshhi$thresh[1]
} else if(FDRthreshlo$thresh[1] > 0){
  FDRthreshold = FDRthreshlo$thresh[1]
} else {
  FDRthreshold = NA
}


FDRthreshhi2 <- phyllos_gwas %>%
  as_tibble() %>%
  filter(between(.data$FDR_adj, 0.01001, 1)) %>%  # 1% FDR threshold
  summarise(thresh = max(.data$log10p))
FDRthreshlo2 <- phyllos_gwas %>%
  as_tibble() %>%
  filter(between(.data$FDR_adj, 0, 0.00999)) %>%  # 1% FDR threshold
  summarise(thresh = min(.data$log10p))
if(FDRthreshhi2$thresh[1] > 0 & FDRthreshlo2$thresh[1] > 0){
  FDRthreshold2 = (FDRthreshhi$thresh[1] + FDRthreshlo$thresh[1])/2
} else if(FDRthreshhi$thresh[1] > 0){
  FDRthreshold2 = FDRthreshhi$thresh[1]
} else if(FDRthreshlo$thresh[1] > 0){
  FDRthreshold2 = FDRthreshlo$thresh[1]
} else {
  FDRthreshold2 = NA
}

#boxplot(-pvals)
#trim for plotting. Trust me you don't want all these points. 
gwas_trim<-as.data.frame(phyllos_gwas[which(pvals<(-3)),])
  gwas_trim$POS<-POS[which(pvals<(-3))]
  gwas_trim$CHR<-CHRN[which(pvals<(-3))]

high_thresh<-.05/nrow(phyllos_gwas)

gwas_trim$CHR2<-paste(substr(gwas_trim$CHR,5,6),sep="")

a<-ggplot(gwas_trim, aes(x=POS, y=(-pvals)))+ 
  geom_point(aes(col=CHR), size=.6)+
  geom_hline(aes(yintercept=-log10(high_thresh)), lty="dashed")+
  geom_hline(aes(yintercept=FDRthreshold))+
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=0, ymax=3.3, fill=gwas_trim$CHR), size=1)+
  scale_color_manual(values = rep(c("black","grey40"),9), guide="none")+
  scale_fill_manual(values = rep(c("black","grey40"),9), guide="none")+
  facet_grid(.~CHR2, scales = "free", space="free")+
  labs(x="Position", y="-log10 pval", 
       title="Microbiome Structure GWA")+
  theme_classic()+
  coord_cartesian(ylim = c(0,11))+
  theme(axis.text.x = element_blank())
a
ggsave("../figs/phyllos_gwas_mar15.png", height=4, width=12)
ggsave("../figs/phyllos_gwas_mar15.pdf", height=4, width=12)
ggsave("../figs/phyllos_gwas_mar15.eps", height=4, width=12)

phyllos_peaks<-gwas_trim[gwas_trim$pvals<(-top10),]
twon<-phyllos_gwas[phyllos_gwas$CHR=="Chr02N",]
twon_peak<-twon[twon$POS>60880000 & twon$POS<60990000,]

b<-ggplot(twon_peak, aes(x=POS, y=(-pvals)))+ 
  geom_point(size=.6)+
  geom_hline(aes(yintercept=-log10(high_thresh)), lty="dashed")+
  geom_hline(aes(yintercept=FDRthreshold))+
  labs(x="Position", y="-log10 pval")+
  annotate("segment", x=60920998, xend=60922038, y=10.5, yend=10.5, col="red", size=1,arrow = arrow(length = unit(0.1, "cm")), lineend="round")+
  annotate("text", x=60936000,y=12.4,label="Pavir.2NG521906", col="red", size=3.5)+
  annotate("text", x=60946000,y=11.7,label="Pavir.2NG521912", col="red", size=3.5)+
  annotate("text", x=60956000,y=11,label="Pavir.2NG521915", col="red", size=3.5)+
  annotate("segment", x=60943081, xend=60946763, y=10.5, yend=10.5, col="red", size=1,arrow = arrow(length = unit(0.1, "cm")), lineend="round")+
  annotate("segment", x=60932143, xend=60936150, y=10.5, yend=10.5, col="red", size=1,arrow = arrow(length = unit(0.1, "cm")), lineend="round")+
  theme_classic()
b
ggsave("../figs/genes_2N_feb28.png", height=4, width=3)
ggsave("../figs/genes_2N_feb28.pdf", height=4, width=3)
ggsave("../figs/genes_2N_feb28.eps", height=4, width=3)


write_csv(phyllos_peaks, "phyllos_gwas_peaks_feb1.csv")  

#checks####
#Chr02N_60928981
#check the frequency of your top SNP.
test15<-snp_subset(x = snp, ind.col = 2345500:2347000)
test15_snp<-snp_attach(test15)
str(test15_snp)
ids<-test15_snp$map$marker.ID
tail(ids)
mat2<-test15_snp$genotypes[]
colnames(mat2)<-test15_snp$map$marker.ID
rownames(mat2)<-test15_snp$fam$sample.ID
write_csv(as.data.frame(mat2),file = "mat2_feb10.csv")

t15<-snp_subset(x = test15_snp, ind.row =test15_snp$fam$sample.ID %in% phenos_gwas2$PLANT_ID )
t15_snp<-snp_attach(t15)
str(t15_snp)

#calculate MAF
mafs<-snp_MAF(t15_snp$genotypes)
names(mafs)<-t15_snp$map$marker.ID

mafs["Chr02N_60928981"]

mat1<-t15_snp$genotypes[]
colnames(mat1)<-t15_snp$map$marker.ID
rownames(mat1)<-t15_snp$fam$sample.ID

write_csv(as.data.frame(mat1),file = "mat1_feb10.csv")

#check for missing data
miss<-mat1[,"Chr02N_60928981"]

table(miss)
gens<-test15_snp$genotypes[]
miss

pheno_gwas4<-phenos_gwas3[phenos_gwas3$PLANT_ID %in% names(miss),]
pheno_gwas5<-pheno_gwas4 %>% group_by(PLANT_ID) %>%
  summarise(mdsum=mean(MDS2))

comp_gens<-data.frame(geno=miss,
                      pheno=pheno_gwas5$mdsum,
                      PLANT_ID=pheno_gwas5$PLANT_ID)
PHEN_META_TRIM_update <- read_csv("~/Desktop/PHEN_META_TRIM_update.csv")

comp_gens$pops<-PHEN_META_TRIM_update$SUBPOP_SNP[match(comp_gens$PLANT_ID, PHEN_META_TRIM_update$PLANT_ID)]

comp_key<-data.frame(pops=unique(comp_gens$pops),
                     Populations=c("Gulf", "Midwest", "Atlantic_admix", "Gulf_admix","Atlantic", "Gulf", "Gulf_admix", "Midwest_admix"))

comp_gens$Subpopulations<-comp_key$Populations[match(comp_gens$pops, comp_key$pops)]

ggplot(comp_gens, aes(x=geno, y=pheno))+
  geom_jitter(aes(col=Subpopulations), width=.2, size=2)+
  scale_color_manual(values=c("cornflowerblue", "blue", "indianred3", "salmon", "purple4", "darkorchid"))+
  labs(x="Genotype", y="NMDS2")+
  theme_classic()
ggsave("../figs/genotype_split_outlier_apr27.png", height = 5, width = 5)
ggsave("../figs/genotype_split_outlier_apr27.pdf", height = 5, width = 5)

DOE <- read_excel("~/Downloads/DOE_GWAS_Master Plant List (2).xlsx")
PVDIV <- read_excel("~/Downloads/PVDIV_Master Metadata File_5-8-2019.xlsx")

#Add in RNA output and plot####
load("~/Desktop/GWAS/rdata/transcript_feb15.RData")

foc_long2$Ecotype<-c("Lowland \n (Resistant)","Upland \n (Susceptible)")[match(foc_long2$eco, c("Lowland", "Upland"))]
foc_long2$Site<-as.character(foc_long2$Site)
foc_long2$Site[foc_long2$Site=="Pickle"]<-"Austin"
foc_long2$Site<-factor(foc_long2$Site, levels = c("KBS", "Columbia", "Austin"))

c<-ggplot(foc_long2, aes(x=geno, y=value, fill=Ecotype))+
  geom_boxplot(fill=NA)+
  geom_jitter(width=.2,  aes(col=Ecotype), size=.6)+
  scale_fill_manual(values=c("firebrick", "skyblue"), name="Ecotype")+
  scale_color_manual(values=c("firebrick", "skyblue"), name="Ecotype")+
  labs(x="Genotype", y="Transcript Counts")+
  facet_grid(id~Site, scales="free")+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle = 90))
c
ggsave("../figs/rna_boxes2_mar15.png", height = 2.25,width=3.5)
ggsave("../figs/rna_boxes2_mar15.pdf", height = 2.25,width=3.5)
ggsave("../figs/rna_boxes2_mar15.eps", height = 2.25,width=3.5)

patch1<-(a/(b+c))


layout<-c(
     patchwork::area(t = 1, l = 1, b = 3, r = 12),
     patchwork::area(t = 4, l = 1, b = 7, r = 12),
     patchwork::area(t = 4, l = 6, b = 7, r = 12)
   )

patch1 + plot_annotation(tag_levels = 'a')+   plot_layout(design=layout)

ggsave("../figs/phyllos_gwas_comb_mar15.png", height = 5,width = 7.5)
ggsave("../figs/phyllos_gwas_comb_mar15.pdf", height = 5,width = 7.5)
ggsave("../figs/phyllos_gwas_comb_mar15.eps", height = 5,width = 7.5)
 