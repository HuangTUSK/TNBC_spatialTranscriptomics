
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

###Differential expresion analysis
dds <- DESeq(dds)
ntd <- normTransform(dds)

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

mtx <- assay(ntd)[sel,]
df <- df[order(df$X),]
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

pheatmap(mtx, cluster_rows=T, show_rownames=T,
         show_colnames=F, cluster_cols=T, fontsize_row = 8,
         annotation_col=df,
         annotation_colors = ann_col)

out <- pheatmap(mtx, cluster_rows=T, show_rownames=T,
                show_colnames=F, cluster_cols=T, fontsize_row = 8)

clust <- data.frame(colnames(mtx[,out$tree_col[["order"]]]))

write.table(clust, "cluster_CD8cell.txt", sep = "\t")

####CD8Tcell sugrroup 
Cd8 <- read.table("cluster_CD8cell.txt", sep = "\t", header = T)
Cd8$segment
mergedDF$SegmentDisplayName

all_merge <- merge(mergedDF, Cd8, by.x="SegmentDisplayName", by.y="segment")
colnames(all_merge)

sub_merge <- all_merge[,c("SegmentDisplayName", "type","Patient","Tissue",
                          "MetType","Survival","TumorType","Label","M2","CD8Tcell")]

write.table(sub_merge,"ImmmuSubtype.txt", sep = "\t", quote = F, row.names = F)


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

#pdf("SpacialType_TNBCsubtype_stackedplot.pdf", width = 6, height = 3)
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
#dev.off() 

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
gne <- read.table("gene.txt", sep = "\t", header = T)
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



#################immune cell compare 
subty <- read.table("Tumor_subtype_immunecellFraction.txt", sep = "\t", header = T)
head(subty)

TNBC <- read.table("samplelist.txt", sep = "\t" , header = T)
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




  
  

