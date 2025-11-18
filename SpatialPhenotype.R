
library(tidyverse)
library(standR)
library(DESeq2)
##meta data
TNBC <- read.table("samplelist.txt", sep = "\t" , header = T)
head(TNBC)
TNBC$DSPID
TNBC$AOI <- gsub("\\..*", "", TNBC$AOI)
TNBC$Label <- gsub("\\_.*", "", TNBC$Label)
TNBC$Patient

TNBC$merge <- paste(TNBC$TNBCtype, TNBC$X, sep="_")
sub11<- TNBC[TNBC$type=="Tumor",]
sub11[sub11$TNBCtype=="M",]

all <- read.table("./TNBCmetaData.txt", sep = "\t", header = T) 
head(all)
all$AOISurfaceArea <- gsub("\\..*", "", all$AOISurfaceArea)

sel <- all[,c("SegmentDisplayName","AOISurfaceArea")]

mergedDF <- merge(sel, TNBC, by.x = "AOISurfaceArea", by.y = "AOI")

rownames(mergedDF) <- mergedDF$SegmentDisplayName

###only stroma
#mergedDF <- mergedDF[mergedDF$Tissue=="LN",]
#mergedDF <- mergedDF[mergedDF$type=="Stroma",]
mergedDF$TumorType <- gsub("-","_", mergedDF$TumorType)
mergedDF$TumorType <- gsub("\\+","_P", mergedDF$TumorType)
###count 
countFile <- read_tsv("allCount.txt") %>% as.data.frame()
head(countFile)
countFile <- countFile[countFile$TargetName!="NegProbe-WTX", ]
rownames(countFile) <- countFile$TargetName
countFile$TargetName <- NULL
###order the rowname and colname
coldata <- countFile[,mergedDF$SegmentDisplayName]

##check 
all(rownames(mergedDF) %in% colnames(coldata))
all(rownames(mergedDF) == colnames(coldata))


##input 
dds <- DESeqDataSetFromMatrix(countData = coldata,
                              colData = mergedDF,
                              design = ~ type) ##type ##Tissue ##M2
dds
##filter 
smallestGroupSize <- 8##for tumor 20
keep <- rowSums(counts(dds) >= 5) >= smallestGroupSize
dds <- dds[keep,]

#factor levels
dds$type <- factor(dds$type, levels = c("Tumor","Stroma"))
#dds$TumorType <- factor(dds$TumorType, levels = c("tnbc","er_pr_her_P"))
dds$M2 <- factor(dds$M2, levels = c("high","low"))
#dds$Tissue <- factor(dds$Tissue, levels = c("LN","breast"))
###Differential expresion analysis
dds <- DESeq(dds)
#res <- results(dds, contrast=c("TumorType","tnbc","er_pr_her_P"))
#res <- results(dds, contrast=c("Tissue","LN","breast"))
res <- results(dds, contrast=c("M2","high","low"))
res <- results(dds, contrast=c("type","Tumor","Stroma"))
#res

resOrdered <- res[order(res$pvalue),]
resSig <- subset(resOrdered, padj < 0.05)

#write.csv(as.data.frame(resSig), file="Tumor_vs_Stroma_results.csv")

###heatmap
ntd <- normTransform(dds)
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
###top 200 DEGs
sel <- rownames(resOrdered[1:50,])
#sel <- rownames(resSig)
#
#sel <- read.table("PAM50GeneList.txt", sep = "\t", header = F)
sel <- read.table("Tcell_marker.txt", sep = "\t", header = F)
sel$V1 <- unique(sel$V1)
#sel <- sel[sel$V1 %in% c("ARG1","ARG2","CD163","MRC1"),]
sel <- c("CD163","PDCD1","CD14","IDO1",
         "ARG1","FOXP3","GZMB")
sel <- c("HAVCR2","CTLA4","CD274","IDO1","IL10","PDCD1","VTCN1","CD276")
sel <- c("CD3E","CD68",'CD8A',"CD8B","MS4A1","NCAM1") ##CD8marker
sel <- c("CD3E","CD68","CD163","MRC1","ARG1",'CD8A',"CD8B","MS4A1","NCAM1") 
###annotation
df <- as.data.frame(colData(dds)[,c("MetType","Label","P53","X","type")])  ###,"TNBCtype.x")])
df$P53Type <- ifelse(df$P53 == "+", "+", "-")
df$P53 <- NULL
#df <- df[order(df$type),]
df <- df[order(df$MetType),]
df$MetType <- sub("\\..*","", df$MetType)
df$MetType <- gsub('[0-9]+', '', df$MetType)
df$MetType <- gsub(' ', '', df$MetType)
#df$Label <-  substr(df$Label, 1, nchar(df$Label)-1)
#df$condition <- paste(df$Patient,df$MetType, df$type, sep = "") 
#df <- df[order(df$condition),]
#df$condition

##order annotation
#sel %in% rownames(assay(ntd))
#sel$V1 %in% rownames(assay(ntd))
#sel_fil <- sel[sel$V1 %in% rownames(assay(ntd)),]
#sel$V1
#sel_fil <- sel[sel$V1!="CD95L",] ##for Tcell
#sel_fil <- sel[!sel$V1 %in% c("AL2","CLB54"),] ##for HIF

#sel_uni <- unique(sel_fil)
#mtx <- assay(ntd)[sel_uni,]
mtx <- assay(ntd)[sel,]
#mtx <- assay(ntd)[sel,]
df <- df[order(df$X),]
#df <- df[df$type=="Tumor",]
mtx <- mtx[,rownames(df)]

##change annotation color 
ann_col = list(
  Patient = c(Slide1_Patient1="black", Slide1_Patient10="grey",
              Slide1_Patient2="red",Slide1_Patient3="yellow",
              Slide1_Patient4="green",Slide1_Patient5="blue",
              Slide1_Patient6="pink",Slide1_Patient7="sienna",
              Slide1_Patient8="magenta",Slide1_Patient9="orange",
              Slide2_Patient10="royalblue",Slide2_Patient13="olivedrab",
              Slide2_Patient19="khaki",Slide2_Patient20="palegreen2",
              Slide2_Patient22="orchid",
              Slide2_Patient28="tan",Slide2_Patient5="cadetblue2",
              Slide2_Patient8="cornsilk",Slide2_Patient9="hotpink2"))

ann_col = list(
  Label = c(P1="black", P10="grey",
            P2="red",P3="yellow",
            P4="green",P5="blue",
            P6="pink",P7="sienna",
            P8="magenta",P9="orange",
            P11="royalblue",P12="olivedrab",
            P13="khaki",P14="palegreen2",
            P15="orchid",P16="tan",
            P17="cadetblue2",P18="cornsilk",
            P19="hotpink2",P20="gold2"))


df$Label
pheatmap(mtx, cluster_rows=T, show_rownames=T,
         show_colnames=F, cluster_cols=T, fontsize_row = 8,
         annotation_col=df,
         annotation_colors = ann_col)

out <- pheatmap(mtx, cluster_rows=T, show_rownames=T,
                show_colnames=F, cluster_cols=T, fontsize_row = 8)
clust <- data.frame(colnames(mtx[,out$tree_col[["order"]]]))

write.table(clust, "cluster_CD8cell.txt", sep = "\t")

####CD8Tcell sugrroup 
Cd8 <- read.table("CD8TcellGroup.txt", sep = "\t", header = T)
Cd8$segment
mergedDF$SegmentDisplayName

all_merge <- merge(mergedDF, Cd8, by.x="SegmentDisplayName", by.y="segment")
colnames(all_merge)

sub_merge <- all_merge[,c("SegmentDisplayName", "type","Patient","Tissue",
                          "MetType","Survival","TumorType","Label","M2","CD8Tcell")]

write.table(sub_merge,"ImmmuSubtype.txt", sep = "\t", quote = F, row.names = F)

spatialtype <- read.table("ImmmuSubtype.txt", sep = "\t", header = T)
head(spatialtype)

####survivalfor all102 samples, p <=0.03
spatialtype$MetType <- sub("\\..*","", spatialtype$MetType)
spatialtype$MetType <- gsub('[0-9]+', '', spatialtype$MetType)
spatialtype$MetType <- gsub(' ', '', spatialtype$MetType)

spatialtype$months <- as.numeric(gsub("[^0-9]", "", spatialtype$Survival))
spatialtype$vital_status <- sub("\\(.*", "", spatialtype$Survival)
spatialtype$vital_status <- gsub(" ","", spatialtype$vital_status)
spatialtype$status[spatialtype$vital_status=="ALIVE"]=0
spatialtype$status[spatialtype$vital_status=="DEAD"]=1
sfit= survfit(Surv(months,status)~MetType,data=spatialtype)
Pvalue=surv_pvalue(sfit)$pval
Pvalue

#####exclude nonmeta
spatialtype$TNBC_spacial <- paste(spatialtype$MetType,spatialtype$SpatialType,sep="_")
spatialtype_meta <- spatialtype[!spatialtype$MetType=="Nonmetas",]
spatialtype_meta <- spatialtype[!spatialtype$SpatialType=="Excluded",]

spatialtype_meta <- spatialtype[!spatialtype$TNBC_spacial=="Nonmetas_Excluded",]
sfit= survfit(Surv(months,status)~TNBCtype,data=spatialtype_meta)  
sfit$
Pvalue=surv_pvalue(sfit)$pval
Pvalue

##plot 
ggsurvplot(sfit,pval=TRUE,
           size=1,
           palette = c("#00a087", "#3c5488", "#f39b7f"),
           break.x.by = 1000,
           title = "",#ggtheme=custom_theme(),
           legend.title = "IFNGR1",
           #legend.labs = c("High", "Low"),
           pval.coord = c(3200, 0.85),
           pval.size = 6,
           xlab="Time / Days",
           font.x = c(20,face="bold"),
           font.y = c(20, face = "bold"),
           font.tickslab = c(18),
           font.legend = c(18),
           font.title=c(22, "bold"))


par(mfrow = c(1, 1))
plot(main = "", sfit,
     col =c("#00a087", "#3c5488", "#f39b7f","wheat4","hotpink3","slategray4"),lty = 1,lwd = 3,
     xlab = "Time (months)",ylab = "Probability of Survival")
legend("bottomleft",
       fill = c("#00a087", "#3c5488", "#f39b7f","wheat4","hotpink3","slategray4"),
       legend = c("LN", "Nonmetas","Primary"),bty = "n") ##"Excluded", "Ignored","Inflamed"

#####################cox hazards risk
spatialtype$SpatialType
res.cox <- coxph(Surv(months, status) ~ SpatialType, data = spatialtype)
res.cox
summary(res.cox)
ggforest(res.cox, data=spatialtype)
### Create nice regression output
model1 <- tbl_regression(res.cox, intercept = TRUE)
as_gt(model1) %>% 
  gt::tab_header("Cox proportional hazards regression model") %>% 
  gt::tab_options(table.align='left')

#####score 
setwd("/Users/emmahuang/Downloads/Bedi_project")

#BiocManager::install("msigdbr")
library(GSVA)
library(msigdbr)
TPM <- read.table("TPMnormalized.txt", sep = "\t", header = T)
head(TPM)
rownames(TPM) <- TPM$X
TPM$X <- NULL
colnames(TPM)
colnames(TPM) <- gsub("\\.","", colnames(TPM))

TPM[1:5,1:5]


gslst <- read.table("gdlist.txt",sep = "\t", header = T)
hallmarks_list <- split(
  gslst$genename, # The genes we want split into pathways
  gslst$geneset # The pathways made as the higher levels of the list
)

gsvaPar <- gsvaParam(data.matrix(TPM), hallmarks_list)
gsva.es <- gsva(gsvaPar, verbose=FALSE)
dim(gsva.es)
res <- data.frame(gsva.es[1:4, 1:112])
tres <- data.frame(t(res))
tres$segment <- rownames(tres)
rownames(tres) <- NULL
spatialtype$SegmentDisplayName <- gsub('\\|', '', spatialtype$SegmentDisplayName)
spatialtype$SegmentDisplayName <- gsub(' ', '', spatialtype$SegmentDisplayName)
score <- merge(spatialtype, tres, by.x="SegmentDisplayName", by.y="segment")

colnames(score)

resDF<-data.frame()
for (column in colnames(score)[c(15:18)]) {
  resDF<-rbind (
    resDF,
    data.frame(Signature=column,
               pvalue=wilcox.test(score[score$SpatialType=="Inflamed",][,column], 
                                  score[score$SpatialType=="Excluded",][,column])$p.value,
               Highmean= mean(score[score$SpatialType=="Inflamed",][,column], na.rm=TRUE),
               Lownmean= mean(score[score$SpatialType=="Excluded",][,column], na.rm=TRUE)
    )
  )
}

resDF


##plot
library(ggsignif)
library(ggplot2)
library(tidyverse)
library(ggpubr)

df_new <- tidyr::gather(score,key='Cell',value='value',15:18)

p1<-ggplot(data = df_new, aes(x=Cell, y=value)) +
  geom_boxplot(aes(fill=SpatialType), outlier.size = 0.1) +
  # ylim(-1,7) +
  labs(fill = "") +
  # ylab("Expression (log10)") +
  ylab("Immune cell score")+
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
  theme(axis.text.y = element_text(size=15, color = "black"),
        axis.text.x = element_text(size=15, color = "black", angle = 45, hjust=1),
        axis.title.y = element_text(size=15, color = "black"),# axis.title = element_text(size=22),
        axis.title.x = element_blank(),plot.title = element_text(size = 27, hjust = 0.5),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        #legend.title = element_text(size=13),
        #legend.position = c(0.9,0.7)
  ) #c(0.9,0.8)
p1 
p1 + stat_compare_means(aes(group = M2),method="wilcox", label = "p.format",
                        label.y = 1, size =6, hide.ns = TRUE) + ggtitle("")
?stat_compare_means


######TNBC subtype # Get the stacked barplot
spatialtype <- read.table("ImmmuSubtype.txt", sep = "\t", header = T)
head(spatialtype)

spatialtype$MetType <- sub("\\..*","", spatialtype$MetType)
spatialtype$MetType <- gsub('[0-9]+', '', spatialtype$MetType)
spatialtype$MetType <- gsub(' ', '', spatialtype$MetType)

subty <- read.table("Tumor_subtype_immunecellFraction.txt", sep = "\t", header = T)
head(subty)
colnames(subty)
subty <- subty[,c(1,7,8,10)]
rownames(subty) <- subty$SegmentDisplayName
colnames(subty)
subty$SegmentDisplayName <- gsub("\\|","",subty$SegmentDisplayName)
subty$SegmentDisplayName <- gsub(" ","",subty$SegmentDisplayName)

subtype_merge <- merge(subty, spatialtype, by.x="SegmentDisplayName", by.y="SegmentDisplayName")
subtype_merge$type
frequency_table <- table(subtype_merge$TNBCtype, subtype_merge$X)
print(frequency_table)

##data
type <- c("BL1","BL2","IM","LAR","M","MSL","UNS")
Excluded <- c(8,4,3,7,0,4,11)
Ignored <- c(0,0,0,3,0,4,3)
Inflamed <- c(0,0,5,2,1,0,1)
mydata <- data.frame(type, Excluded, Ignored,Inflamed)
library(tidyverse)
mydata <- mydata %>% pivot_longer(cols = !type, names_to = "SpatialType", values_to = "Frequencies") 

pdf("../Figures/SpacialType_TNBCsubtype_stackedplot.pdf", width = 6, height = 3)
  ggplot(mydata, aes(fill = SpatialType, x = type, y = Frequencies)) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = c("tan","cadetblue2","#3c5488")) +
  xlab("") +  
  theme_bw() +
  theme(axis.line = element_line(color = 'black')) +
  theme(axis.text.y = element_text(size=14, color = "black"),
        axis.title.y = element_text(size=14, color = "black"),
        axis.text.x = element_text(size=12, color = "black"),
          # axis.text.x = element_text(size=15, color = "black"), axis.title = element_text(size=22),
        axis.title.x = element_blank(),plot.title = element_text(size = 14, hjust = 0.5),
        legend.text = element_text(size=14),
        legend.title = element_blank())
dev.off() 

par(mfrow=c(1, 1), mar=c(5, 5, 4, 4))                                
barplot(frequency_table, 
        col=colors()[c(23,89,12,50,7,3,40)] , 
        border="white", 
        space=0.04, 
        font.axis=2, 
        xlab="",
        legend.text = TRUE,
        args.legend = list(x = "topright", bty = "n", inset=c(-0.2, 0)))


#########workflow plot alluvial plot
TNBC <- read.table("TNBCalllist.txt", sep = "\t" , header = T)
head(TNBC)

lst <- read_tsv("./TNBCsample.txt") %>% as.data.frame()
head(lst)
nrow(lst)
lst$area <- gsub("\\..*", "", lst$area)
lst <- lst[lst$Sample_ID %in% TNBC$DSPID,]

all <- read.table("./TNBCmetaData.txt", sep = "\t", header = T) 
head(all)
all$SegmentDisplayName <- gsub("\\|","",all$SegmentDisplayName)
all$SegmentDisplayName <- gsub(" ","",all$SegmentDisplayName)

spatialtype <- read.table("ImmmuSubtype.txt", sep = "\t", header = T)
head(spatialtype)

spatialtype$MetType <- sub("\\..*","", spatialtype$MetType)
spatialtype$MetType <- gsub('[0-9]+', '', spatialtype$MetType)
spatialtype$MetType <- gsub(' ', '', spatialtype$MetType)

spatialtype$SegmentDisplayName <- gsub("\\|","",spatialtype$SegmentDisplayName)
spatialtype$SegmentDisplayName <- gsub(" ","",spatialtype$SegmentDisplayName)

subty <- read.table("Tumor_subtype_immunecellFraction.txt", sep = "\t", header = T)
head(subty)
colnames(subty)
subty <- subty[,c(1,7,8,10)]
rownames(subty) <- subty$SegmentDisplayName
colnames(subty)
subty$SegmentDisplayName <- gsub("\\|","",subty$SegmentDisplayName)
subty$SegmentDisplayName <- gsub(" ","",subty$SegmentDisplayName)

subtype_merge <- merge(subty, spatialtype, by.x="SegmentDisplayName", by.y="SegmentDisplayName")
subtype_merge$SegmentDisplayName

all_merge <- merge(all, subtype_merge,by.x="SegmentDisplayName",by.y="SegmentDisplayName")
all_merge$AOISurfaceArea <- gsub("\\..*", "", all_merge$AOISurfaceArea)

sampleAnnoFile <- all_merge[all_merge$AOISurfaceArea %in% lst$area,]

head(sampleAnnoFile)
nrow(sampleAnnoFile)
sampleAnnoFile$SegmentDisplayName

###count 
countFile <- read_tsv("allCount.txt") %>% as.data.frame()
head(countFile)
colnames(countFile) <- gsub("\\|","",colnames(countFile))
colnames(countFile) <- gsub(" ","",colnames(countFile))

countFile2 <- countFile[, sampleAnnoFile$SegmentDisplayName]
countFile2$TargetName <- countFile$TargetName

###featureAnnoFile
gne <- read.table("/Users/emmahuang/Downloads/Bedi_project/gene.txt", sep = "\t", header = T)
head(gne)
gsub("\\,.*","",gne$location)
gne$start <- as.numeric(sub(".*:(.*)-.*", "\\1", gne$location))
gne$end <- as.numeric(sub(".*-", "", gne$location))

gne$GeneLength <- gne$end - gne$start
nrow(gne)
featureAnnoFile <- read_tsv("GeneMeta.txt") %>% as.data.frame()
head(featureAnnoFile)
nrow(featureAnnoFile)
head(sampleAnnoFile)
head(countFile2)

##
library(ggalluvial)
spe <- readGeoMx(countFile2, sampleAnnoFile, featureAnnoFile,colnames.as.rownames = c("TargetName", "SegmentDisplayName", "TargetName"))
spe@colData$CD8Tcell

standR::plotSampleInfo(spe, column2plot = c("SpatialType","Label","TNBCtype","MetType"))
?plotSampleInfo


###############################DEseq2 
##meta data
TNBC <- read.table("./Deseq2/samplelist.txt", sep = "\t" , header = T)
head(TNBC)
TNBC$DSPID
TNBC$AOI <- gsub("\\..*", "", TNBC$AOI)
TNBC$Label <- gsub("\\_.*", "", TNBC$Label)

table(TNBC$Tissue,TNBC$type, TNBC$X)

all <- read.table("./TNBCmetaData.txt", sep = "\t", header = T) 
head(all)
all$AOISurfaceArea <- gsub("\\..*", "", all$AOISurfaceArea)

sel <- all[,c("SegmentDisplayName","AOISurfaceArea")]

mergedDF <- merge(sel, TNBC, by.x = "AOISurfaceArea", by.y = "AOI")

rownames(mergedDF) <- mergedDF$SegmentDisplayName

###only stroma
mergedDF <- mergedDF[mergedDF$X %in% c("Excluded","Inflamed"),]
mergedDF <- mergedDF[mergedDF$Tissue=="breast",]
mergedDF <- mergedDF[mergedDF$type=="Tumor",]
mergedDF$TumorType <- gsub("-","_", mergedDF$TumorType)
mergedDF$TumorType <- gsub("\\+","_P", mergedDF$TumorType)
###count 
countFile <- read_tsv("allCount.txt") %>% as.data.frame()
head(countFile)
countFile <- countFile[countFile$TargetName!="NegProbe-WTX", ]
rownames(countFile) <- countFile$TargetName
countFile$TargetName <- NULL
###order the rowname and colname
coldata <- countFile[,mergedDF$SegmentDisplayName]

##check 
all(rownames(mergedDF) %in% colnames(coldata))
all(rownames(mergedDF) == colnames(coldata))


##input 
dds <- DESeqDataSetFromMatrix(countData = coldata,
                              colData = mergedDF,
                              design = ~ M2) ##type ##Tissue ##M2
dds
##filter 
smallestGroupSize <- 5##for tumor 20
keep <- rowSums(counts(dds) >= 5) >= smallestGroupSize
dds <- dds[keep,]

#factor levels
dds$X <- factor(dds$X, levels = c("Excluded","Inflamed")) #"Excluded","Inflamed"
#dds$TumorType <- factor(dds$TumorType, levels = c("tnbc","er_pr_her_P"))
dds$M2 <- factor(dds$M2, levels = c("high","low"))
###Differential expresion analysis
dds <- DESeq(dds)
res <- results(dds, contrast=c("X","Excluded","Inflamed"))
#res <- results(dds, contrast=c("Tissue","LN","breast"))
res <- results(dds, contrast=c("M2","high","low"))
#res <- results(dds, contrast=c("type","Tumor","Stroma"))
#res

resOrdered <- res[order(res$pvalue),]
resSig <- subset(resOrdered, padj < 0.05)

write.csv(as.data.frame(resSig), file="Excluded_vs_Inflamed_DEG.csv")

###heatmap
ntd <- normTransform(dds)
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
###top 200 DEGs
sel <- rownames(resOrdered[1:50,])
#sel <- rownames(resSig)
#
#sel <- read.table("PAM50GeneList.txt", sep = "\t", header = F)
sel <- read.table("M2marker.txt", sep = "\t", header = F)
sel$V1 <- unique(sel$V1)
#sel <- sel[sel$V1 %in% c("ARG1","ARG2","CD163","MRC1"),]
sel2 <- c(sel$V1,"CTLA4","PDCD1","CD274","LAG3","HAVCR2","TIGIT", "BTLA","VTCN1")
sel <- c("HAVCR2","CTLA4","CD274","PDCD1","IDO1","IL10","CD276","VTCN1")
sel <- c("CD3E","CD68",'CD8A',"CD8B","MS4A1","NCAM1") ##CD8marker
sel <- c("TMSB15A","CALML5","SDC1","GPRC5A","CXADR","PBX1","IGFBP5",
         "ENTPD3","TUFT1","GREM1","COL10A1","THBS2","PERP",#"CEACAM6",
         "TCN1","SCGB2A2","SCGB1A1","AKR1C2","MMP2","CPE","ASPN",
         "TPSAB1","COL1A1","SPOCK1","FAP","CAMK2N1","SPON1","COL5A1",#"WARS",
         "CXCL13","CCL5","GZMB","TRBC1","LCK","CORO1A","SIRPG","PLAC8",
         "PVRIG","CCL18","IL7R","IL2RG","NKG7","IGHG1")##paper markers
###annotation
df <- as.data.frame(colData(dds)[,c("MetType","Label","P53","X","M2","type")])  ###,"TNBCtype.x")])
df$P53Type <- ifelse(df$P53 == "+", "+", "-")
df$P53 <- NULL
#df <- df[order(df$type),]
df$MetType <- sub("\\..*","", df$MetType)
df$MetType <- gsub('[0-9]+', '', df$MetType)
df$MetType <- gsub(' ', '', df$MetType)
#df$Label <-  substr(df$Label, 1, nchar(df$Label)-1)
#df$condition <- paste(df$Patient,df$MetType, df$type, sep = "") 
#df <- df[order(df$condition),]
#df$condition

##order annotation
sel %in% rownames(assay(ntd))
#sel$V1 %in% rownames(assay(ntd))
#sel_fil <- sel[sel$V1 %in% rownames(assay(ntd)),]
#sel$V1
#sel_fil <- sel[sel$V1!="CD95L",] ##for Tcell
#sel_fil <- sel[!sel$V1 %in% c("AL2","CLB54"),] ##for HIF

#sel_uni <- unique(sel_fil)
#mtx <- assay(ntd)[sel$V1,]
#mtx <- assay(ntd)[sel2,]
##filter 
df <- df[df$type=="Tumor",]
mtx <- assay(ntd)[sel,]
df$type <- NULL
df <- df[order(df$X),]
#df <- df[df$type=="Tumor",]
mtx <- mtx[,rownames(df)]

##change annotation color 
ann_col = list(
  Patient = c(Slide1_Patient1="black", Slide1_Patient10="grey",
              Slide1_Patient2="red",Slide1_Patient3="yellow",
              Slide1_Patient4="green",Slide1_Patient5="blue",
              Slide1_Patient6="pink",Slide1_Patient7="sienna",
              Slide1_Patient8="magenta",Slide1_Patient9="orange",
              Slide2_Patient10="royalblue",Slide2_Patient13="olivedrab",
              Slide2_Patient19="khaki",Slide2_Patient20="palegreen2",
              Slide2_Patient22="orchid",
              Slide2_Patient28="tan",Slide2_Patient5="cadetblue2",
              Slide2_Patient8="cornsilk",Slide2_Patient9="hotpink2"))

ann_col = list(
  Label = c(P1="black", P10="grey",
            P2="red",P3="yellow",
            P4="green",P5="blue",
            P6="pink",P7="sienna",
            P8="magenta",P9="orange",
            P11="royalblue",P12="olivedrab",
            P13="khaki",P14="palegreen2",
            P15="orchid",P16="tan",
            P17="cadetblue2",P18="cornsilk",
            P19="hotpink2",P20="gold2"),
  M2 = c(high="hotpink2", low="gold2"),
  X = c(Excluded="grey62",Ignored="tan",Inflamed="cadetblue2"),
  MetType = c(LN="royalblue",Nonmetas="orchid", Primary="palegreen2"))


df$Label
pheatmap(mtx, cluster_rows=F, show_rownames=T,
         show_colnames=F, cluster_cols=F, fontsize_row = 8,
         annotation_col=df,
         annotation_colors = ann_col)


####box plot 
tmtx <- data.frame(t(mtx))
tmtx$sample <- rownames(tmtx)

tmtx_merge<- merge(tmtx, mergedDF, by.x="sample", by.y="SegmentDisplayName")

######
##tumor
tmtx_merge$Spatial_TNBC <- paste(tmtx_merge$TNBCtype,tmtx_merge$X, sep ="_")
colnames(tmtx_merge)

resDF<-data.frame()
for (column in colnames(tmtx_merge)[c(2:9)]) {
  resDF<-rbind (
    resDF,
    data.frame(Signature=column,
              # tst <- aov(tmtx_merge[,column] ~ tmtx_merge[,61]),
               pvalue=summary(aov(tmtx_merge[,column] ~ tmtx_merge[,32]))[[1]][["Pr(>F)"]][1],
             #  pvalue=wilcox.test(tmtx_merge[tmtx_merge$X=="Inflamed",][,column], 
             #                     tmtx_merge[tmtx_merge$X=="Ignored",][,column])$p.value,
               Excluded_mean= mean(tmtx_merge[tmtx_merge$X=="Excluded",][,column], na.rm=TRUE),
             Ignored_mean= mean(tmtx_merge[tmtx_merge$X=="Ignored",][,column], na.rm=TRUE),
               Inflamed_mean= mean(tmtx_merge[tmtx_merge$X=="Inflamed",][,column], na.rm=TRUE)
    )
  )
}

resDF
resDF[resDF$pvalue<0.05,]
#write.table(resDF, "Anova_spatialtypeDEG.txt", sep = "\t", quote = F, row.names = F)
###plot the markers 
resDF2 <- resDF[!resDF$Signature %in% c("COL1A1","IGHG1","CAMK2N1"),]

heatmapDF <- resDF2[,c(3,4,5)]

rownames(heatmapDF) <- resDF2$Signature
heatmapDF_log <- log2(heatmapDF)
pheatmap(heatmapDF_log, cluster_rows=F, 
         show_colnames=F, cluster_cols=F, fontsize_row = 8)

##plot
library(ggsignif)
library(ggplot2)
library(tidyverse)
library(ggpubr)

df_new <- tidyr::gather(tmtx_merge,key='Gene',value='exp',2:9)

p1<-ggplot(data = df_new, aes(x=Gene, y=exp)) +
  geom_boxplot(aes(fill=TNBCtype), outlier.size = 0.1) +
  # ylim(-1,7) +
  labs(fill = "") +
  ylab("Expression (log10)") +
  #ylab("Immune cell score")+
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
    values = c("#e64b35", "#4dbbd5","#00a087", "#3c5488", "#f39b7f","wheat4","hotpink3","slategray4",
                                "royalblue","olivedrab",
                                "tan","cadetblue2","yellow3")) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black')) +
  theme(axis.text.y = element_text(size=15, color = "black"),
        axis.text.x = element_text(size=15, color = "black"),
        # axis.text.x = element_text(size=15, color = "black"), axis.title = element_text(size=22),
        axis.title.x = element_blank(),plot.title = element_text(size = 27, hjust = 0.5),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        #legend.title = element_text(size=13),
        #legend.position = c(0.9,0.7)
  ) #c(0.9,0.8)
p1 

#pdf("./Figures/CheckpointGene_Spatial_TNBC_boxplot.pdf", width = 10,height = 5)
pdf("./Figures/CheckpointGene_TNBC_boxplot.pdf", width = 10,height = 5)

p1 + stat_compare_means(aes(group = TNBCtype),method="anova", label = "p.signif",
                        label.y = 7, size =6, hide.ns = TRUE) + ggtitle("")

dev.off()

tumor <- df_new[df_new$type=="Tumor",]
tumor$MetType <- sub("\\..*","", tumor$MetType)
tumor$MetType <- gsub('[0-9]+', '', tumor$MetType)
tumor$MetType <- gsub(' ', '', tumor$MetType)
tumor <- tumor[tumor$MetType %in% c("Primary","Nonmetas"),]

####plot 
p1<-ggplot(data = tumor, aes(x=Gene, y=exp)) +
  geom_boxplot(aes(fill=MetType), outlier.size = 0.1) +
  # ylim(-1,7) +
  labs(fill = "") +
  ylab("Expression (log10)") +
  #ylab("Immune cell score")+
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
    values = c("#e64b35", "#4dbbd5","#00a087", "#3c5488", "#f39b7f","wheat4","hotpink3","slategray4",
               "royalblue","olivedrab",
               "tan","cadetblue2","yellow3")) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black')) +
  theme(axis.text.y = element_text(size=15, color = "black"),
        axis.text.x = element_text(size=15, color = "black"),
        # axis.text.x = element_text(size=15, color = "black"), axis.title = element_text(size=22),
        axis.title.x = element_blank(),plot.title = element_text(size = 27, hjust = 0.5),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        #legend.title = element_text(size=13),
        #legend.position = c(0.9,0.7)
  ) #c(0.9,0.8)
p1 

#pdf("./Figures/CheckpointGene_Spatial_TNBC_boxplot.pdf", width = 10,height = 5)
pdf("./Figures/CheckpointGene_MetNonMet_boxplot.pdf", width = 10,height = 5)

p1 + stat_compare_means(aes(group = MetType),method="wilcox", label = "p.signif",
                        label.y = 7, size =6, hide.ns = TRUE) + ggtitle("")

dev.off()

###PCA
library("RColorBrewer")
vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$type, vsd$MetType, sep="-")
colnames(sampleDistMatrix) <- NULL

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

plotPCA(vsd, intgroup=c("MetType"))
colData(dds)
pcaData <- plotPCA(vsd, intgroup=c("Label","X", "Tissue"), returnData=TRUE)

pcaData$MetType <- sub("\\..*","", pcaData$MetType)
pcaData$MetType <- gsub('[0-9]+', '', pcaData$MetType)
pcaData$MetType <- gsub(' ', '', pcaData$MetType)
pcaData$Met_type <- paste(pcaData$MetType, pcaData$type,sep = "_")

percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Label, shape = X)) +
  geom_point(size=3,aes(shape=factor(X))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(size=4) +
  scale_shape_manual(values=1:15) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  scale_color_manual(values = c("P1"="black", "P10"="grey",
                                "P2"="red","P3"="yellow",
                                "P4"="green","P5"="blue",
                                "P6"="pink","P7"="sienna",
                                "P8"="magenta","P9"="orange",
                                "P11"="royalblue","P12"="olivedrab",
                                "P13"="khaki","P14"="palegreen2",
                                "P15"="orchid","P16"="tan",
                                "P17"="cadetblue2","P18"="cornsilk",
                                "P19"="hotpink2","P20"="gold2"))


#####################################GSEA
#prepare the list for fgsea
res$SYMBOL <- rownames(res)
res<- data.frame(res)
res2 <- res %>% 
  dplyr::select(SYMBOL, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(SYMBOL) %>% 
  summarize(stat=mean(stat))
res2

#devtools::install_github("ctlab/fgsea")
library(fgsea)
ranks <- deframe(res2)
head(ranks, 20)

#########################################
#look at KEGG pathways:
pathways.kegg <- gmtPathways("../msigdb_v2024.1.Hs_files_to_download_locally/msigdb_v2024.1.Hs_GMTs/c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt")

fgseaResT <- fgsea(pathways=pathways.kegg, stats=ranks)
head(fgseaResT[order(pval), ])

sum(fgseaResT[, pval < 0.05])

#fwrite(fgseaResT, file="CR_vs_ICR_fgseaRes_kegg.txt", sep="\t", sep2=c("", " ", ""))

#plot the GSEA results

#gene_list<- read.table("CR_vs_ICR_fgseaRes_kegg.txt",sep="\t",header=TRUE)
gene_list <- data.frame(fgseaResT)
gene_list<- gene_list[with(gene_list,order(NES)),]

top1<- gene_list[1:50,]
top1 

gene_list<- gene_list[with(gene_list,order(-NES)),]

top2<- gene_list[1:50,]
top2 
top50 <- rbind(top1, top2)
head(top50,3)
select<- names(top50) %in% c("pathway","pval","padj","NES" ) 
df_for_lollipop_plot <- top50[select]
#change pathway to lowercase
df_for_lollipop_plot <- df_for_lollipop_plot %>% mutate(pathway = tolower(pathway))
##delete hallmark characters 
df_for_lollipop_plot <- df_for_lollipop_plot %>% mutate(pathway = str_sub(pathway,6,-1))
names(df_for_lollipop_plot)<-c( "Pathway","pval","padj","NES")
#re-order plot
df_for_lollipop_plot<-df_for_lollipop_plot[order(df_for_lollipop_plot$NES), ]
print(df_for_lollipop_plot)
df_for_lollipop_plot$Pathway <- factor(df_for_lollipop_plot$Pathway, levels = df_for_lollipop_plot$Pathway)
df_for_lollipop_plot$pval<-as.numeric(df_for_lollipop_plot$pval)
#filter by p-value and ratio
df_for_lollipop_plot <- subset(df_for_lollipop_plot, pval < 0.05)
# df_for_lollipop_plot<-df_for_lollipop_plot[(df_for_lollipop_plot$pval < 0.1)]
df_for_lollipop_plot$Pathway
#######filter some pathway
#fil_lst <- c("asthma", "allograft_rejection","graft_versus_host_disease",
#             "type_i_diabetes_mellitus", "autoimmune_thyroid_disease","asthma", "viral_myocarditis",
#             "non_small_cell_lung_cancer", "pancreatic_cancer","acute_myeloid_leukemia"). ## for UNS

fil_lst <- c("small_cell_lung_cancer", "prion_diseases","bladder_cancer",
             "glioma", "leishmania_infection","prostate_cancer", "arrhythmogenic_right_ventricular_cardiomyopathy_arvc",
             "chronic_myeloid_leukemia", "pathways_in_cancer","endometrial_cancer",
             "thyroid_cancer","pathogenic_escherichia_coli_infection") ## for LAR

fil_lst <- c("small_cell_lung_cancer", "prion_diseases","bladder_cancer","pathogenic_escherichia_coli_infection",
             "glioma", "thyroid_cancer","non_small_cell_lung_cancer","acute_myeloid_leukemia",
             "pancreatic_cancer","chronic_myeloid_leukemia","type_i_diabetes_mellitus","asthma",
             "leishmania_infection","prostate_cancer", "arrhythmogenic_right_ventricular_cardiomyopathy_arvc",
             "chronic_myeloid_leukemia", "pathways_in_cancer","endometrial_cancer",
             "thyroid_cancer","pathogenic_escherichia_coli_infection") ## for BL

fil_lst <- c("asthma", "allograft_rejection","autoimmune_thyroid_disease",
             "type_i_diabetes_mellitus", "viral_myocarditis","graft_versus_host_disease", 
             "leishmania_infection","parkinsons_disease",
             "huntingtons_disease", "small_cell_lung_cancer","non_small_cell_lung_cancer",
             "prion_diseases","chronic_myeloid_leukemia","pancreatic_cancer",
             "endometrial_cancer","alzheimers_disease") ## for IM

fil_lst <- c("asthma", "allograft_rejection","autoimmune_thyroid_disease",
             "type_i_diabetes_mellitus", "viral_myocarditis","graft_versus_host_disease", 
             "leishmania_infection","parkinsons_disease","glioma","pathways_in_cancer",
             "huntingtons_disease", "small_cell_lung_cancer","non_small_cell_lung_cancer",
             "prion_diseases","chronic_myeloid_leukemia","pancreatic_cancer",
             "endometrial_cancer","alzheimers_disease","colorectal_cancer",
             "acute_myeloid_leukemia","vibrio_cholerae_infection",
             "terpenoid_backbone_biosynthesis","pathogenic_escherichia_coli_infection",
             "systemic_lupus_erythematosus","renal_cell_carcinoma",
             "vasopressin_regulated_water_reabsorption","long_term_potentiation",
             "progesterone_mediated_oocyte_maturation","maturity_onset_diabetes_of_the_young",
             "drug_metabolism_cytochrome_p450","terpenoid_backbone_biosynthesis","proteasome",
             "ribosome","lysosome","spliceosome","antigen_processing_and_presentation",
             "prostate_cancer","thyroid_cancer") ## for IM

df_for_lollipop_plot <- df_for_lollipop_plot[!df_for_lollipop_plot$Pathway %in% fil_lst,]

df_for_lollipop_plot[[ "log_pvalue" ]]<-as.numeric(-log10(df_for_lollipop_plot$pval))
df_for_lollipop_plot[[ "color_scheme" ]]<-ifelse(df_for_lollipop_plot$NES < 0, "red", "blue")
df_for_lollipop_plot <- df_for_lollipop_plot[order(df_for_lollipop_plot$NES), ]
my_cols<-as.factor(df_for_lollipop_plot[[ "color_scheme" ]])
maxsize<-max(df_for_lollipop_plot$log_p)

p1<-ggplot(df_for_lollipop_plot, aes(x=Pathway, y=NES, label=NES)) + 
  geom_point(stat='identity', aes(col=my_cols, size=log_pvalue))+ 
  scale_colour_manual("NES", values=c("red", "blue"), labels=c("Up", "Down"))+
  labs(y = "NES", x = "KEGG Pathway") +
  scale_size_continuous(range = c(2, 8.5)) +  # <-  this will have the effect of pulling the dot sizes towards one another. may need to comment this line for future data sets.
  geom_segment(aes(y = 0, x = Pathway, yend = NES, xend = Pathway, col=my_cols))+
  ylim(c(-3, 3)) + coord_flip() +
  theme(axis.text.x = element_text(color="black", size=12, angle=0),
        axis.text.y = element_text(color="black", size=11, angle=0)) +
  theme( axis.line = element_line(colour = "black", size = 1, linetype = "solid"))+
  theme(axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"))+ 
  theme(legend.title = element_text(colour="black", size=14, face="bold"))+
  theme(legend.text = element_text(colour="black", size=14, face="bold"))+
  geom_hline(yintercept=0, linetype="dashed", color = "gray")+
  theme(
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    plot.background = element_blank()) 

p1

#################immune cell compare 
subty <- read.table("Tumor_subtype_immunecellFraction.txt", sep = "\t", header = T)
head(subty)

TNBC <- read.table("./Deseq2/samplelist.txt", sep = "\t" , header = T)
head(TNBC)
TNBC$DSPID
TNBC$AOI <- gsub("\\..*", "", TNBC$AOI)
TNBC$Label <- gsub("\\_.*", "", TNBC$Label)

all <- read.table("./TNBCmetaData.txt", sep = "\t", header = T) 
head(all)
all$AOISurfaceArea <- gsub("\\..*", "", all$AOISurfaceArea)

sel <- all[,c("SegmentDisplayName","AOISurfaceArea")]

mergedDF <- merge(sel, TNBC, by.x = "AOISurfaceArea", by.y = "AOI")
colnames(mergedDF)
mergedDF2 <- mergedDF[,c(2,19,20,21,22,23)]

Immune_merge <- merge(mergedDF2,subty,by.x="SegmentDisplayName",by.y="SegmentDisplayName")
colnames(Immune_merge)
Immune_merge$MetType <- sub("\\..*","", Immune_merge$MetType)
Immune_merge$MetType <- gsub('[0-9]+', '', Immune_merge$MetType)
Immune_merge$MetType <- gsub(' ', '', Immune_merge$MetType)

Immune_merge <- Immune_merge[!Immune_merge$MetType=="LN",] ##"Primary"  "LN"  "Nonmetas" 
resDF<-data.frame()
for (column in colnames(Immune_merge)[c(30:51)]) {
  resDF<-rbind (
    resDF,
    data.frame(Signature=column,
               # tst <- aov(tmtx_merge[,column] ~ tmtx_merge[,61]),
               pvalue=summary(aov(Immune_merge[,column] ~ Immune_merge[,4]))[[1]][["Pr(>F)"]][1],
               #  pvalue=wilcox.test(tmtx_merge[tmtx_merge$X=="Inflamed",][,column], 
               #                     tmtx_merge[tmtx_merge$X=="Ignored",][,column])$p.value,
               Excluded_mean= mean(Immune_merge[Immune_merge$X=="Excluded",][,column], na.rm=TRUE),
               Ignored_mean= mean(Immune_merge[Immune_merge$X=="Ignored",][,column], na.rm=TRUE),
               Inflamed_mean= mean(Immune_merge[Immune_merge$X=="Inflamed",][,column], na.rm=TRUE)
    )
  )
}

resDF


df_new <- tidyr::gather(Immune_merge,key='Gene',value='exp',30:51)

p1<-ggplot(data = df_new, aes(x=Gene, y=exp)) +
  geom_boxplot(aes(fill=X), outlier.size = 0.1) +
  # ylim(-1,7) +
  labs(fill = "") +
  # ylab("Expression (log10)") +
  ylab("Immune cell score")+
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
  theme(axis.text.y = element_text(size=15, color = "black"),
        axis.text.x = element_text(size=10, color = "black", angle = 45, hjust=1),
        # axis.text.x = element_text(size=15, color = "black"), axis.title = element_text(size=22),
        axis.title.x = element_blank(),plot.title = element_text(size = 27, hjust = 0.5),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        #legend.title = element_text(size=13),
        #legend.position = c(0.9,0.7)
  ) #c(0.9,0.8)
p1 
p1 + stat_compare_means(aes(group = X),method="anova", label = "p.signif",
                        label.y = 0.32, size =6, hide.ns = TRUE) + ggtitle("")




  
  

