
library(SpatialDecon)
library(tidyverse)
library(dplyr)

##input prepare 
norm = read.table("TPMnormalized.txt", sep = "\t", header = T)
head(norm)
rownames(norm) <- norm$X
norm$X <- NULL

##annotation file
TNBC <- read.table("TNBCalllist.txt", sep = "\t" , header = T)
head(TNBC)

lst <- read_tsv("./TNBCsample.txt") %>% as.data.frame()
head(lst)
nrow(lst)
lst$area <- gsub("\\..*", "", lst$area)
lst <- lst[lst$Sample_ID %in% TNBC$DSPID,]
lst$Organ

tissue <- lst[,c("Organ","area")]

all <- read.table("./TNBCmetaData.txt", sep = "\t", header = T) 
head(all)
all$AOISurfaceArea <- gsub("\\..*", "", all$AOISurfaceArea)
sampleAnnoFile <- all[all$AOISurfaceArea %in% lst$area,]

sampleAnnoFile2 <- merge(sampleAnnoFile, tissue, by.x="AOISurfaceArea", by.y = "area")

head(sampleAnnoFile2)
nrow(sampleAnnoFile2)
colnames(sampleAnnoFile2)

Anno <- sampleAnnoFile2[,c(2,4,5,6,11,29,30,33,34)]

head(Anno)

###
Anno$Label <- paste("ROI",Anno$ROILabel, Anno$Type, sep = "")
Anno$ROI <- paste("ROI", Anno$ROILabel, sep = "")
Anno$AOI.name <- Anno$Type
Anno$x <- Anno$ScanOffsetX
Anno$y <- Anno$ScanOffsetY
Anno$nuclei <- Anno$AOINucleiCount
Anno$Sample_ID <- Anno$SegmentDisplayName
Anno$Silde.name <- Anno$SlideName

Anno <- Anno %>%
  mutate(AOI.annotation=case_when(Type
    =="Stroma" ~ "TME",
    Type =="Tumor"~ "BreastTumor"))
Anno$Tissue <- Anno$Organ
rownames(Anno) <- Anno$SegmentDisplayName
Anno <- Anno %>%
  mutate(istumor=case_when(Type =="Stroma" ~ "FALSE",
                           Type =="Tumor"~ "TRUE"))

head(Anno)

colnames(Anno)
Anno2 <- Anno[,c(16,20,17,11,12,18,13,14,15,19)]

###count 
countFile <- read_tsv("allCount.txt") %>% as.data.frame()
head(countFile)
colnames(countFile)
countFile2 <- countFile[, sampleAnnoFile$SegmentDisplayName]
countFile$TargetName= make.names(countFile$TargetName, unique = T)
rownames(countFile2) <- countFile$TargetName


####
colnames(norm)  
colnames(countFile2) 
rownames(Anno2)
colnames(norm) <- gsub("\\."," ",colnames(norm))
colnames(norm) <- gsub("  "," |",colnames(norm) )

### use the NegProbe to estimate per-observation background
per.observation.neg = norm[grep('^NegProbe', rownames(norm)),]

per.observation.mean.neg = colMeans(per.observation.neg)

# and define a background matrix in which each column (observation) is the
# appropriate value of per-observation background:???
bg = sweep(norm * 0, 2, per.observation.mean.neg, "+")
dim(bg)

##For human, also have other tissue type
data("safeTME")
data("safeTME.matches")
signif(safeTME[seq_len(3), seq_len(3)], 2)

heatmap(sweep(safeTME, 1, apply(safeTME, 1, max), "/"),
        labRow = NA, margins = c(10, 5))

##human breast cancer 
humanbreast <- download_profile_matrix(species = "Human",
                                       age_group = "Adult", 
                                       matrixname = "BreastCancer_Wu") #
dim(mousespleen)
head(cellGroups)
metadata
heatmap(sweep(mousespleen, 1, apply(mousespleen, 1, max), "/"),
        labRow = NA, margins = c(10, 5), cexCol = 0.7)

##
is.matrix(norm)
norm.matrix <- data.matrix(norm)

countFile2 <- countFile2[rownames(norm),]
countFile2 <- na.omit(countFile2)
norm <- norm[rownames(countFile2),]
norm.matrix <- data.matrix(norm)
countFile2.matrix <- data.matrix(countFile2)

res = spatialdecon(norm = norm.matrix, 
                   raw = countFile2.matrix,
                   bg = bg,
                   X = humanbreast,#humanbreast, safeTME
                   align_genes = TRUE)

##Performing basic deconvolution with the spatialdecon function
res$se
write.table(res$beta, "Breast_cancer_cell_tumor.txt", sep = "\t", quote = F, row.names = T)

###“beta”, the matrix of estimated cell abundances.
heatmap(t(res$beta), cexCol = 0.5, cexRow = 0.7, margins = c(10,7))

##advanced setting of spatialdecon
str(safeTME.matches)
Anno2$istumor = (Anno2$AOI.name == "Tumor")

na.omit(norm)

countFile2 <- countFile2[rownames(norm),]
countFile2 <- na.omit(countFile2)
norm <- norm[rownames(countFile2),]
norm.matrix <- data.matrix(norm)
countFile2.matrix <- data.matrix(countFile2)
ncol(countFile2.matrix)
ncol(norm.matrix)
ncol(bg)
restils = spatialdecon(norm = norm.matrix,                     # normalized data
                       raw = countFile2.matrix,                       # raw data, used to down-weight low-count observations 
                       bg = bg,                         # expected background counts for every data point in norm
                       X = safeTME,                     # safeTME matrix, used by default
                    #   cellmerges = humanbreast.matches,   # safeTME.matches object, used by default
                       cell_counts = Anno2$nuclei,      # nuclei counts, used to estimate total cells
                    #   is_pure_tumor = Anno2$istumor,   # identities of the Tumor segments/observations
                       n_tumor_clusters = 10)   

restils$cell.counts
heatmap(t(res$beta), cexCol = 0.5, cexRow = 0.7, margins = c(10,7))
write.table(res$beta, "Breast_cancer_cellabundancescore_tumor.txt", sep = "\t", quote = F, row.names = T)


###plot
data("cellcols")
o = hclust(dist(t(restils$cell.counts$cell.counts)))$order
layout(mat = (matrix(c(1, 2), 1)), widths = c(34, 34))
TIL_barplot(restils$cell.counts$cell.counts[, o], draw_legend = TRUE, 
            cex.names = 0.5)

TIL_barplot(mat = res$beta,draw_legend=TRUE)

?layout
?TIL_barplot

# or the proportions of cells:
temp = replace(restils$prop_of_nontumor, is.na(restils$prop_of_nontumor), 0)
o = hclust(dist(temp[restils$AOI.name == "TME",]))$order
TIL_barplot(t(restils$prop_of_nontumor[restils$AOI.name == "TME",])[, o], 
            draw_legend = TRUE, cex.names = 0.75)


####following analysis

####result from spatialDecon
setwd("/Users/emmahuang/Downloads/Bedi_project")

TNBC <- read.table("./Deseq2/samplelist.txt", sep = "\t" , header = T)
head(TNBC)
TNBC$DSPID
TNBC$AOI <- gsub("\\..*", "", TNBC$AOI)

all <- read.table("./TNBCmetaData.txt", sep = "\t", header = T) 
head(all)
all$AOISurfaceArea <- gsub("\\..*", "", all$AOISurfaceArea)

sel <- all[,c("SegmentDisplayName","AOISurfaceArea")]

mergedDF <- merge(sel, TNBC, by.x = "AOISurfaceArea", by.y = "AOI")

mergedDF$SegmentDisplayName <- gsub(" \\| ", "_", mergedDF$SegmentDisplayName)
mergedDF$SegmentDisplayName <- gsub(" ", "", mergedDF$SegmentDisplayName)

###spatialDecon breast cell component
####result from spatialDecon
setwd("/Users/emmahuang/Downloads/Bedi_project")

TNBC <- read.table("./Deseq2/samplelist.txt", sep = "\t" , header = T)
head(TNBC)
TNBC$DSPID
TNBC$AOI <- gsub("\\..*", "", TNBC$AOI)

all <- read.table("./TNBCmetaData.txt", sep = "\t", header = T) 
head(all)
all$AOISurfaceArea <- gsub("\\..*", "", all$AOISurfaceArea)

sel <- all[,c("SegmentDisplayName","AOISurfaceArea")]

mergedDF <- merge(sel, TNBC, by.x = "AOISurfaceArea", by.y = "AOI")

mergedDF$SegmentDisplayName <- gsub(" \\| ", "_", mergedDF$SegmentDisplayName)
mergedDF$SegmentDisplayName <- gsub(" ", "", mergedDF$SegmentDisplayName)

###spatialDecon breast cell component
TME <- read.table("Breast_cancer_cell_tumor.txt", sep = "\t", header = T)
head(TME)

rowSums(TME)
###remove the cell had sum is 0 
notlst <- c("Myeloid_c5_Macrophage_3_SIGLEC1", "Myeloid_c11_cDC2_CD1C", 
            "Myeloid_c9_Macrophage_2_CXCL10","Myeloid_c8_Monocyte_2_S100A9", 
            "Myeloid_c1_LAM1_FABP5", "Myeloid_c2_LAM2_APOE", 
            "Myeloid_c12_Monocyte_1_IL1B","Myeloid_c10_Macrophage_1_EGR1", 
            "T_cells_c10_NKT_cells_FCGR3A", "T_cells_c9_NK_cells_AREG",
            "T_cells_c3_CD4+_Tfh_CXCL13", "T_cells_c2_CD4+_T-regs_FOXP3",
            "T_cells_c8_CD8+_LAG3", "T_cells_c7_CD8+_IFNG",
            "T_cells_c5_CD8+_GZMK","T_cells_c4_CD8+_ZFP36",
            "Endothelial.CXCL12","Endothelial.ACKR1")

TME <- subset(TME, !rownames(TME) %in% notlst)

t(TME)
tTME <- data.frame(t(TME))

ncol(tTME)
tTME$Sum <- rowSums(tTME)

###transform to fraction 
for (i in 1:31) {
  tTME[[i]] <- tTME[[i]] / tTME[[32]]
}

tTME$Sum <- NULL

rownames(tTME) <- gsub("\\."," ",rownames(tTME))
rownames(tTME) <- gsub("  ","_",rownames(tTME))
rownames(tTME) <- gsub(" ","",rownames(tTME))
tTME$sampleID <- rownames(tTME)
rownames(tTME) <- NULL

mergedCancerCell <- merge(mergedDF, tTME, by.x = "SegmentDisplayName", by.y = "sampleID")

####compare 
colnames(mergedCancerCell)

resDF<-data.frame()
for (column in colnames(mergedCancerCell)[c(18:48)]) {
  resDF<-rbind (
    resDF,
    data.frame(Signature=column,
               pvalue=wilcox.test(mergedCancerCell[mergedCancerCell$Tissue=="breast",][,column], 
                                  mergedCancerCell[mergedCancerCell$Tissue=="LN",][,column])$p.value,
               Tumormean= mean(mergedCancerCell[mergedCancerCell$Tissue=="breast",][,column], na.rm=TRUE),
               Stromamean= mean(mergedCancerCell[mergedCancerCell$Tissue=="LN",][,column], na.rm=TRUE)
    )
  )
  
  
}

resDF

##plot
library(ggsignif)
library(ggplot2)
library(tidyverse)
library(ggpubr)

df_new <- tidyr::gather(mergedCancerCell,key='cellType',value='Fraction',18:48,-1)
df_new$type

p1<-ggplot(data = df_new, aes(x=cellType, y=Fraction)) +
  geom_boxplot(aes(fill=type), outlier.size = 0.1) +
   ylim(0,1.5) +
  labs(fill = "") +
  # ylab("Expression (log10)") +
  ylab("Cancer cells components")+
  theme_bw() +
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  scale_fill_manual(
    #labels = c(expression(IFNGR1^High),expression(IFNGR1^Low)),
 #   values = c("#e64b35", "#4dbbd5","#00a087", "#3c5488", "#f39b7f","wheat4","hotpink3","slategray4")) +
    values = c("#3c5488", "#f39b7f"))+
    #draws x and y axis line
  theme(axis.line = element_line(color = 'black')) +
  theme(axis.text.y = element_text(size=14, color = "black"),
        axis.title.y = element_text(size=14, color = "black"),
        axis.text.x = element_text(size=12, color = "black", angle = 45, hjust=1),
        # axis.text.x = element_text(size=15, color = "black"), axis.title = element_text(size=22),
        axis.title.x = element_blank(),plot.title = element_text(size = 14, hjust = 0.5),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        #legend.title = element_text(size=13),
        legend.position = c(0.9,0.7)
  ) #c(0.9,0.8)
p1 
p1 + stat_compare_means(aes(group = type),method="wilcox.test", label = "p.signif",
                        label.y = 1.2, size =6, hide.ns = TRUE) + ggtitle("")


###stack bar plot 
tumor_set <- df_new[df_new$type=="Stroma",]
ggplot(df_new, aes(fill=cellType, y=Fraction, x=DSPID)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("") +
  xlab("") +
  ylab("Cell Components") +
  theme_bw() +
  theme(axis.text.y = element_text(size=12, color = "black"),
        axis.title = element_text(size=14, color = "black"),
        axis.text.x = element_text(size=5, color = "black", angle = 45, hjust = 1, vjust = 1),
        #axis.text.x = element_blank(),
        #axis.text = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=10, color = "black"),
        panel.spacing.x = unit(1, "mm"),
        strip.background.x = element_blank(),
        strip.text.x = element_blank()) +
  facet_grid(.~Patient, scales = "free_x") +
  scale_fill_manual(values = c("hotpink","hotpink3",
                               "khaki","khaki1","khaki2","khaki3","khaki4",
                               "skyblue","skyblue1","skyblue2","skyblue3","skyblue4",
                               "sienna1","sienna3",
                               "#00a081","#00a089",
                               "#3c5488", 
                               "#f39b7f",
                               "darkseagreen","darkseagreen1","darkseagreen3","darkseagreen4",
                               "thistle",
                               "purple",
                               "azure1","azure2","azure3",
                               "slategray1","slategray2","slategray3","slategray4"
  ))



######subtype compare 
TNBCtype <- read.table("./TNBCtype/TNBC_sutype_result.csv", sep = ",", header = T)
head(TNBCtype)

TNBCtype$X <- gsub("\\.+", "_", gsub("\\.", "", TNBCtype$X))
mergedDF$SegmentDisplayName <- gsub("_","", mergedDF$SegmentDisplayName)
mergedsubtype <- merge(mergedDF, TNBCtype, by.x = "SegmentDisplayName", by.y = "X")

tTME$sampleID <- gsub("_","", tTME$sampleID)
Immune_subtype <- merge(mergedsubtype, tTME, by.x = "SegmentDisplayName", by.y = "sampleID")

colnames(Immune_subtype)

df_new <- tidyr::gather(Immune_subtype,key='cellType',value='Fraction',21:51,-1)
df_new$subtype

tumorset <- df_new[df_new$type=="Tumor",]

p1<-ggplot(data = tumorset, aes(x=cellType, y=Fraction)) +
  geom_boxplot(aes(fill=subtype), outlier.size = 0.1) +
  ylim(0,1) +
  labs(fill = "") +
  # ylab("Expression (log10)") +
  ylab("Cell components")+
  theme_bw() +
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  scale_fill_manual(
    #labels = c(expression(IFNGR1^High),expression(IFNGR1^Low)),
    values = c("#e64b35", "#4dbbd5","#00a087", "#3c5488", "#f39b7f","wheat4","hotpink3","slategray4")) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black')) +
  theme(axis.text.y = element_text(size=14, color = "black"),
        axis.title.y = element_text(size=14, color = "black"),
        axis.text.x = element_text(size=12, color = "black", angle = 45, hjust=1),
        # axis.text.x = element_text(size=15, color = "black"), axis.title = element_text(size=22),
        axis.title.x = element_blank(),plot.title = element_text(size = 27, hjust = 0.5),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        #legend.title = element_text(size=13),
        #legend.position = c(0.9,0.7)
  ) #c(0.9,0.8)
p1 
p1 + stat_compare_means(aes(group = subtype),method="anova", label = "p.signif",
                        label.y = 0.9, size =6, hide.ns = TRUE) + ggtitle("")

?stat_compare_means

colnames(Immune_subtype)

##
tumorset <- Immune_subtype[Immune_subtype$type=="Tumor",]
colnames(tumorset)

resDF<-data.frame()
for (column in colnames(tumorset)[c(21:51)]) {
  resDF<-rbind (
    resDF,
    data.frame(Signature=column,
               pvalue=wilcox.test(tumorset[tumorset$Tissue=="breast",][,column], 
                            tumorset[tumorset$Tissue=="LN",][,column])$p.value,
               Tumormean= mean(tumorset[tumorset$Tissue=="breast",][,column], na.rm=TRUE),
               Stromamean= mean(tumorset[tumorset$Tissue=="LN",][,column], na.rm=TRUE)
    )
  )
  
  
}

resDF

tumorset <- df_new[df_new$type=="Tumor",]

p1<-ggplot(data = tumorset, aes(x=cellType, y=Fraction)) +
  geom_boxplot(aes(fill=Tissue), outlier.size = 0.1) +
  ylim(0,1.5) +
  labs(fill = "") +
  # ylab("Expression (log10)") +
  ylab("Cancer cells components")+
  theme_bw() +
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  scale_fill_manual(
    #labels = c(expression(IFNGR1^High),expression(IFNGR1^Low)),
    #   values = c("#e64b35", "#4dbbd5","#00a087", "#3c5488", "#f39b7f","wheat4","hotpink3","slategray4")) +
    values = c("#3c5488", "#f39b7f"))+
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black')) +
  theme(axis.text.y = element_text(size=14, color = "black"),
        axis.title.y = element_text(size=14, color = "black"),
        axis.text.x = element_text(size=12, color = "black", angle = 45, hjust=1),
        # axis.text.x = element_text(size=15, color = "black"), axis.title = element_text(size=22),
        axis.title.x = element_blank(),plot.title = element_text(size = 14, hjust = 0.5),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        #legend.title = element_text(size=13),
        legend.position = c(0.9,0.9)
  ) #c(0.9,0.8)
p1 
p1 + stat_compare_means(aes(group = type),method="wilcox.test", label = "p.signif",
                        label.y = 1.2, size =6, hide.ns = TRUE) + ggtitle("")

TME <- read.table("Breast_cancer_cell_tumor.txt", sep = "\t", header = T)
head(TME)

rowSums(TME)
###remove the cell had sum is 0 
notlst <- c("Myeloid_c5_Macrophage_3_SIGLEC1", "Myeloid_c11_cDC2_CD1C", 
            "Myeloid_c9_Macrophage_2_CXCL10","Myeloid_c8_Monocyte_2_S100A9", 
            "Myeloid_c1_LAM1_FABP5", "Myeloid_c2_LAM2_APOE", 
            "Myeloid_c12_Monocyte_1_IL1B","Myeloid_c10_Macrophage_1_EGR1", 
            "T_cells_c10_NKT_cells_FCGR3A", "T_cells_c9_NK_cells_AREG",
            "T_cells_c3_CD4+_Tfh_CXCL13", "T_cells_c2_CD4+_T-regs_FOXP3",
            "T_cells_c8_CD8+_LAG3", "T_cells_c7_CD8+_IFNG",
            "T_cells_c5_CD8+_GZMK","T_cells_c4_CD8+_ZFP36",
            "Endothelial.CXCL12","Endothelial.ACKR1")

TME <- subset(TME, !rownames(TME) %in% notlst)

t(TME)
tTME <- data.frame(t(TME))

ncol(tTME)
tTME$Sum <- rowSums(tTME)

###transform to fraction 
for (i in 1:31) {
  tTME[[i]] <- tTME[[i]] / tTME[[32]]
}

tTME$Sum <- NULL

rownames(tTME) <- gsub("\\."," ",rownames(tTME))
rownames(tTME) <- gsub("  ","",rownames(tTME))
rownames(tTME) <- gsub(" ","",rownames(tTME))
tTME$sampleID <- rownames(tTME)
rownames(tTME) <- NULL

mergedCancerCell <- merge(mergedDF, tTME, by.x = "SegmentDisplayName", by.y = "sampleID")

####compare 
colnames(mergedCancerCell)

resDF<-data.frame()
for (column in colnames(mergedCancerCell)[c(18:48)]) {
  resDF<-rbind (
    resDF,
    data.frame(Signature=column,
               pvalue=wilcox.test(mergedCancerCell[mergedCancerCell$Tissue=="breast",][,column], 
                                  mergedCancerCell[mergedCancerCell$Tissue=="LN",][,column])$p.value,
               Tumormean= mean(mergedCancerCell[mergedCancerCell$Tissue=="breast",][,column], na.rm=TRUE),
               Stromamean= mean(mergedCancerCell[mergedCancerCell$Tissue=="LN",][,column], na.rm=TRUE)
    )
  )
  
  
}

resDF

##plot
library(ggsignif)
library(ggplot2)
library(tidyverse)
library(ggpubr)

df_new <- tidyr::gather(mergedCancerCell,key='cellType',value='Fraction',18:48,-1)
df_new$type

p1<-ggplot(data = df_new, aes(x=cellType, y=Fraction)) +
  geom_boxplot(aes(fill=type), outlier.size = 0.1) +
  ylim(0,1.5) +
  labs(fill = "") +
  # ylab("Expression (log10)") +
  ylab("Cancer cells components")+
  theme_bw() +
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  scale_fill_manual(
    #labels = c(expression(IFNGR1^High),expression(IFNGR1^Low)),
    #   values = c("#e64b35", "#4dbbd5","#00a087", "#3c5488", "#f39b7f","wheat4","hotpink3","slategray4")) +
    values = c("#3c5488", "#f39b7f"))+
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black')) +
  theme(axis.text.y = element_text(size=14, color = "black"),
        axis.title.y = element_text(size=14, color = "black"),
        axis.text.x = element_text(size=12, color = "black", angle = 45, hjust=1),
        # axis.text.x = element_text(size=15, color = "black"), axis.title = element_text(size=22),
        axis.title.x = element_blank(),plot.title = element_text(size = 14, hjust = 0.5),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        #legend.title = element_text(size=13),
        legend.position = c(0.9,0.7)
  ) #c(0.9,0.8)
p1 
p1 + stat_compare_means(aes(group = type),method="wilcox.test", label = "p.signif",
                        label.y = 1.2, size =6, hide.ns = TRUE) + ggtitle("")


###stack bar plot 
tumor_set <- df_new[df_new$type=="Stroma",]
ggplot(df_new, aes(fill=cellType, y=Fraction, x=DSPID)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("") +
  xlab("") +
  ylab("Cell Components") +
  theme_bw() +
  theme(axis.text.y = element_text(size=12, color = "black"),
        axis.title = element_text(size=14, color = "black"),
        axis.text.x = element_text(size=5, color = "black", angle = 45, hjust = 1, vjust = 1),
        #axis.text.x = element_blank(),
        #axis.text = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=10, color = "black"),
        panel.spacing.x = unit(1, "mm"),
        strip.background.x = element_blank(),
        strip.text.x = element_blank()) +
  facet_grid(.~Patient, scales = "free_x") +
  scale_fill_manual(values = c("hotpink","hotpink3",
                               "khaki","khaki1","khaki2","khaki3","khaki4",
                               "skyblue","skyblue1","skyblue2","skyblue3","skyblue4",
                               "sienna1","sienna3",
                               "#00a081","#00a089",
                               "#3c5488", 
                               "#f39b7f",
                               "darkseagreen","darkseagreen1","darkseagreen3","darkseagreen4",
                               "thistle",
                               "purple",
                               "azure1","azure2","azure3",
                               "slategray1","slategray2","slategray3","slategray4"
  ))



######subtype compare 
TNBCtype <- read.table("./TNBCtype/TNBC_sutype_result.csv", sep = ",", header = T)
head(TNBCtype)

TNBCtype$X <- gsub("\\.+", "_", gsub("\\.", "", TNBCtype$X))
mergedDF$SegmentDisplayName <- gsub("_","", mergedDF$SegmentDisplayName)
mergedsubtype <- merge(mergedDF, TNBCtype, by.x = "SegmentDisplayName", by.y = "X")

tTME$sampleID <- gsub("_","", tTME$sampleID)
Immune_subtype <- merge(mergedsubtype, tTME, by.x = "SegmentDisplayName", by.y = "sampleID")

colnames(Immune_subtype)

df_new <- tidyr::gather(Immune_subtype,key='cellType',value='Fraction',21:51,-1)
df_new$subtype

tumorset <- df_new[df_new$type=="Tumor",]

p1<-ggplot(data = tumorset, aes(x=cellType, y=Fraction)) +
  geom_boxplot(aes(fill=subtype), outlier.size = 0.1) +
  ylim(0,1.5) +
  labs(fill = "") +
  # ylab("Expression (log10)") +
  ylab("Cell components")+
  theme_bw() +
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  scale_fill_manual(
    #labels = c(expression(IFNGR1^High),expression(IFNGR1^Low)),
    values = c("#e64b35", "#4dbbd5","#00a087", "#3c5488", "#f39b7f","wheat4","hotpink3","slategray4")) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black')) +
  theme(axis.text.y = element_text(size=14, color = "black"),
        axis.title.y = element_text(size=14, color = "black"),
        axis.text.x = element_text(size=12, color = "black", angle = 45, hjust=1),
        # axis.text.x = element_text(size=15, color = "black"), axis.title = element_text(size=22),
        axis.title.x = element_blank(),plot.title = element_text(size = 27, hjust = 0.5),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        #legend.title = element_text(size=13),
        #legend.position = c(0.9,0.7)
  ) #c(0.9,0.8)
p1 
p1 + stat_compare_means(aes(group = type),method="anova", label = "p.signif",
                        label.y = 1.2, size =12, hide.ns = TRUE) + ggtitle("")

?stat_compare_means

colnames(Immune_subtype)

##
tumorset <- Immune_subtype[Immune_subtype$type=="Tumor",]
colnames(tumorset)

resDF<-data.frame()
for (column in colnames(tumorset)[c(21:51)]) {
  resDF<-rbind (
    resDF,
    data.frame(Signature=column,
               pvalue=wilcox.test(tumorset[tumorset$Tissue=="breast",][,column], 
                                  tumorset[tumorset$Tissue=="LN",][,column])$p.value,
               Tumormean= mean(tumorset[tumorset$Tissue=="breast",][,column], na.rm=TRUE),
               Stromamean= mean(tumorset[tumorset$Tissue=="LN",][,column], na.rm=TRUE)
    )
  )
  
  
}

resDF

tumorset <- df_new[df_new$type=="Tumor",]

p1<-ggplot(data = tumorset, aes(x=cellType, y=Fraction)) +
  geom_boxplot(aes(fill=Tissue), outlier.size = 0.1) +
  ylim(0,1.5) +
  labs(fill = "") +
  # ylab("Expression (log10)") +
  ylab("Cancer cells components")+
  theme_bw() +
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  scale_fill_manual(
    #labels = c(expression(IFNGR1^High),expression(IFNGR1^Low)),
    #   values = c("#e64b35", "#4dbbd5","#00a087", "#3c5488", "#f39b7f","wheat4","hotpink3","slategray4")) +
    values = c("#3c5488", "#f39b7f"))+
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black')) +
  theme(axis.text.y = element_text(size=14, color = "black"),
        axis.title.y = element_text(size=14, color = "black"),
        axis.text.x = element_text(size=12, color = "black", angle = 45, hjust=1),
        # axis.text.x = element_text(size=15, color = "black"), axis.title = element_text(size=22),
        axis.title.x = element_blank(),plot.title = element_text(size = 14, hjust = 0.5),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        #legend.title = element_text(size=13),
        legend.position = c(0.9,0.9)
  ) #c(0.9,0.8)
p1 
p1 + stat_compare_means(aes(group = type),method="wilcox.test", label = "p.signif",
                        label.y = 1.2, size =6, hide.ns = TRUE) + ggtitle("")


########################
TME <- read.table("Breast_cancer_safeTME_tumor.txt", sep = "\t", header = T)
head(TME)

rowSums(TME)

t(TME)
tTME <- data.frame(t(TME))

ncol(tTME)
tTME$Sum <- rowSums(tTME)

###transform to fraction 
for (i in 1:18) {
  tTME[[i]] <- tTME[[i]] / tTME[[19]]
}

tTME$Sum <- NULL

rownames(tTME) <- gsub("\\."," ",rownames(tTME))
rownames(tTME) <- gsub("  ","",rownames(tTME))
rownames(tTME) <- gsub(" ","",rownames(tTME))
tTME$sampleID <- rownames(tTME)
rownames(tTME) <- NULL

mergedCancerCell <- merge(mergedDF, tTME, by.x = "SegmentDisplayName", by.y = "sampleID")

####compare 
colnames(mergedCancerCell)

resDF<-data.frame()
for (column in colnames(mergedCancerCell)[c(24:41)]) {
  resDF<-rbind (
    resDF,
    data.frame(Signature=column,
               pvalue=wilcox.test(mergedCancerCell[mergedCancerCell$type=="Tumor",][,column], 
                                  mergedCancerCell[mergedCancerCell$type=="Stroma",][,column])$p.value,
               Tumormean= mean(mergedCancerCell[mergedCancerCell$type=="Tumor",][,column], na.rm=TRUE),
               Stromamean= mean(mergedCancerCell[mergedCancerCell$type=="Stroma",][,column], na.rm=TRUE)
    )
  )
  
  
}

resDF

##plot
library(ggsignif)
library(ggplot2)
library(tidyverse)
library(ggpubr)

df_new <- tidyr::gather(mergedCancerCell,key='cellType',value='Fraction',24:41,-1)
df_new$type

p1<-ggplot(data = df_new, aes(x=cellType, y=Fraction)) +
  geom_boxplot(aes(fill=type), outlier.size = 0.1) +
  ylim(0,1.5) +
  labs(fill = "") +
  # ylab("Expression (log10)") +
  ylab("Cancer TME")+
  theme_bw() +
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  scale_fill_manual(
    #labels = c(expression(IFNGR1^High),expression(IFNGR1^Low)),
    #   values = c("#e64b35", "#4dbbd5","#00a087", "#3c5488", "#f39b7f","wheat4","hotpink3","slategray4")) +
    values = c("#3c5488", "#f39b7f"))+
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black')) +
  theme(axis.text.y = element_text(size=14, color = "black"),
        axis.title.y = element_text(size=14, color = "black"),
        axis.text.x = element_text(size=12, color = "black", angle = 45, hjust=1),
        # axis.text.x = element_text(size=15, color = "black"), axis.title = element_text(size=22),
        axis.title.x = element_blank(),plot.title = element_text(size = 14, hjust = 0.5),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        #legend.title = element_text(size=13),
        legend.position = c(0.9,0.7)
  ) #c(0.9,0.8)
p1 
p1 + stat_compare_means(aes(group = type),method="wilcox.test", label = "p.signif",
                        label.y = 1.2, size =6, hide.ns = TRUE) + ggtitle("")


###stack bar plot 

mergedCancerCell
mergedCancerCell$figurelabel <- paste(mergedCancerCell$MetType,mergedCancerCell$type,sep = "_")
colnames(mergedCancerCell)

df_new <- tidyr::gather(mergedCancerCell,key='cellType',value='Fraction',24:41,-1)
df_new$type
ggplot(df_new, aes(fill=cellType, y=Fraction, x=figurelabel)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("") +
  xlab("") +
  ylab("Cell Components") +
  theme_bw() +
  theme(axis.text.y = element_text(size=12, color = "black"),
        axis.title = element_text(size=14, color = "black"),
        axis.text.x = element_text(size=5, color = "black", angle = 45, hjust = 1, vjust = 1),
        #axis.text.x = element_blank(),
        #axis.text = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=10, color = "black"),
        panel.spacing.x = unit(1, "mm"),
        strip.background.x = element_blank(),
        strip.text.x = element_blank()) +
  facet_grid(.~Patient, scales = "free_x") +
  scale_fill_manual(values = c("hotpink","hotpink3",
                               "khaki","khaki1","khaki2","khaki3","khaki4",
                               "skyblue","skyblue1","skyblue2","skyblue3","skyblue4",
                               "sienna1","sienna3",
                               "#00a081","#00a089",
                               "#3c5488", 
                               "#f39b7f",
                               "darkseagreen","darkseagreen1","darkseagreen3","darkseagreen4",
                               "thistle",
                               "purple",
                               "azure1","azure2","azure3",
                               "slategray1","slategray2","slategray3","slategray4"
  ))



######subtype compare 
TNBCtype <- read.table("./TNBCtype/TNBC_sutype_result.csv", sep = ",", header = T)
head(TNBCtype)

TNBCtype$X <- gsub("\\.+", "_", gsub("\\.", "", TNBCtype$X))
mergedDF$SegmentDisplayName <- gsub("_","", mergedDF$SegmentDisplayName)
mergedsubtype <- merge(mergedDF, TNBCtype, by.x = "SegmentDisplayName", by.y = "X")

tTME$sampleID <- gsub("_","", tTME$sampleID)
Immune_subtype <- merge(mergedsubtype, tTME, by.x = "SegmentDisplayName", by.y = "sampleID")

colnames(Immune_subtype)

df_new <- tidyr::gather(Immune_subtype,key='cellType',value='Fraction',27:44,-1)
df_new$subtype

tumorset <- df_new[df_new$type=="Tumor",]

p1<-ggplot(data = tumorset, aes(x=cellType, y=Fraction)) +
  geom_boxplot(aes(fill=subtype), outlier.size = 0.1) +
  ylim(0,1) +
  labs(fill = "") +
  # ylab("Expression (log10)") +
  ylab("Cancer TME")+
  theme_bw() +
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  scale_fill_manual(
    #labels = c(expression(IFNGR1^High),expression(IFNGR1^Low)),
    values = c("#e64b35", "#4dbbd5","#00a087", "#3c5488", "#f39b7f","wheat4","hotpink3","slategray4")) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black')) +
  theme(axis.text.y = element_text(size=14, color = "black"),
        axis.title.y = element_text(size=14, color = "black"),
        axis.text.x = element_text(size=12, color = "black", angle = 45, hjust=1),
        # axis.text.x = element_text(size=15, color = "black"), axis.title = element_text(size=22),
        axis.title.x = element_blank(),plot.title = element_text(size = 27, hjust = 0.5),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        #legend.title = element_text(size=13),
        #legend.position = c(0.9,0.7)
  ) #c(0.9,0.8)
p1 
p1 + stat_compare_means(aes(group = subtype),method="anova", label = "p.signif",
                        label.y = 0.9, size =6, hide.ns = T) + ggtitle("")

?stat_compare_means

