
library(tidyverse)
library(standR)
library(DESeq2)
##meta data
TNBC <- read.table("samplelist.txt", sep = "\t" , header = T)
head(TNBC)
TNBC$DSPID
TNBC$AOI <- gsub("\\..*", "", TNBC$AOI)
####summary of clinical factors


######
all <- read.table("./TNBCmetaData.txt", sep = "\t", header = T) 
head(all)
all$AOISurfaceArea <- gsub("\\..*", "", all$AOISurfaceArea)

sel <- all[,c("SegmentDisplayName","AOISurfaceArea")]

mergedDF <- merge(sel, TNBC, by.x = "AOISurfaceArea", by.y = "AOI")

rownames(mergedDF) <- mergedDF$SegmentDisplayName

###only tumor 
#mergedDF <- mergedDF[mergedDF$Tissue=="breast",]
mergedDF <- mergedDF[mergedDF$type=="Tumor",]
mergedDF$MetType <- sub("\\..*","", mergedDF$MetType)
mergedDF$MetType <- gsub('[0-9]+', '', mergedDF$MetType)
mergedDF$MetType <- gsub(' ', '', mergedDF$MetType)
mergedDF <- mergedDF[mergedDF$MetType %in% c("Primary","Nonmetas"),]

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
                              design = ~ MetType) ##type ##Tissue##TumorType ##MetType
dds
##filter 
smallestGroupSize <- 20##for tumor 20
keep <- rowSums(counts(dds) >= 5) >= smallestGroupSize
dds <- dds[keep,]

#factor levels
dds$MetType <- factor(dds$MetType, levels = c("Nonmetas","Primary")) 
#dds$type <- factor(dds$type, levels = c("Tumor","Stroma"))
dds$Tissue <- factor(dds$Tissue, levels = c("LN","breast"))
#dds$TumorType <- factor(dds$TumorType, levels = c("er-pr-her+","tnbc"))
###Differential expresion analysis
dds <- DESeq(dds)
#res <- results(dds, contrast=c("type","Tumor","Stroma"))
#res <- results(dds, contrast=c("TumorType","er-pr-her+","tnbc"))
res <- results(dds, contrast=c("MetType","Nonmetas","Primary"))
res <- results(dds, contrast=c("Tissue","LN","breast"))
res

resOrdered <- res[order(res$pvalue),]
sigset <- subset(resOrdered, pvalue < 0.05)
sigset[sigset$log2FoldChange<0,]

resSig <- subset(resOrdered, padj < 0.05)
resSig
#write.csv(as.data.frame(resSig), file="Tumor_vs_Stroma_results.csv")
###gene plot 

plotCounts(dds, gene="ATRX", intgroup="Tissue")

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
sel <- read.table("Interferon_marker.txt", sep = "\t", header = F)
sel$V1 <- unique(sel$V1)

###annotation
df <- as.data.frame(colData(dds)[,c("MetType","Label","type","TumorType")])  ###,"TNBCtype.x")])
as.data.frame(colData(dds)[,c("MetType","Label","TumorType")])

df$type <- NULL
#df <- df[order(df$type),]
df$MetType <- sub("\\..*","", df$MetType)
df$MetType <- gsub('[0-9]+', '', df$MetType)
df$MetType <- gsub(' ', '', df$MetType)
df$Label <- gsub("\\_.*", "", df$Label)
df <- df[order(df$MetType),]
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
sel <- c("ATRX")

sel_uni <- unique(sel_fil)
mtx <- assay(ntd)[sel_uni,]
mtx <- assay(ntd)[sel,]
#mtx <- assay(ntd)[sel$V1,]
mtx <- mtx[,rownames(df)]

##change annotation color 
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
  MetType = c(LN="#3c5488",Nonmetas="#00a087",Primary="#4dbbd5"))

df$Patient
#pdf("./Figures/NonMetvsPrimary_heatmaptop50.pdf", width = 8,height = 7)
pheatmap(mtx, cluster_rows=T, show_rownames=T,
         show_colnames=FALSE, cluster_cols=F, fontsize_row = 8,
         annotation_col=df,
         annotation_colors = ann_col)

dev.off()

pheatmap(mtx, cluster_rows=T, show_rownames=T,
         show_colnames=FALSE, cluster_cols=F, fontsize_row = 8,
         annotation_col=df,
         annotation_colors = ann_col)
###
library("RColorBrewer")
vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$type, vsd$MetType, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

##only tumor

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

###PCA
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

plotPCA(vsd, intgroup=c("TNBCtype"))

pcaData <- plotPCA(vsd, intgroup=c("MetType","Patient","TumorType","TNBCtype"), returnData=TRUE)

pcaData$MetType <- sub("\\..*","", pcaData$MetType)
pcaData$MetType <- gsub('[0-9]+', '', pcaData$MetType)
pcaData$MetType <- gsub(' ', '', pcaData$MetType)
pcaData$Met_type <- paste(pcaData$MetType, pcaData$TNBCtype,sep = "_")

percentVar <- round(100 * attr(pcaData, "percentVar"))

unique(pcaData$Patient)
pcaData <- pcaData %>%
  mutate(PLabel=case_when
         (pcaData$Patient=="Slide1_Patient10" ~ "P10",
           pcaData$Patient=="Slide1_Patient1" ~ "P1",
           pcaData$Patient=="Slide1_Patient2" ~ "P2",
           pcaData$Patient=="Slide1_Patient3" ~ "P3",
           pcaData$Patient=="Slide1_Patient4" ~ "P4",
           pcaData$Patient=="Slide1_Patient5" ~ "P5",
           pcaData$Patient=="Slide1_Patient6" ~ "P6",
           pcaData$Patient=="Slide1_Patient7" ~ "P7",
           pcaData$Patient=="Slide1_Patient8" ~ "P8",
           pcaData$Patient=="Slide1_Patient9" ~ "P9",
           pcaData$Patient=="Slide2_Patient5" ~ "P11",
           pcaData$Patient=="Slide2_Patient8" ~ "P12",
           pcaData$Patient=="Slide2_Patient9" ~ "P13",
           pcaData$Patient=="Slide2_Patient10" ~ "P14",
           pcaData$Patient=="Slide2_Patient13" ~ "P15",
           pcaData$Patient=="Slide2_Patient19" ~ "P16",
           pcaData$Patient=="Slide2_Patient20" ~ "P17",
           pcaData$Patient=="Slide2_Patient22" ~ "P18",
           pcaData$Patient=="Slide2_Patient28" ~ "P19"
         ))

pcaData$PLabel

pdf("./Figures/Deseq2_PCA.pdf",width = 12, height = 10)

ggplot(pcaData, aes(PC1, PC2, color=PLabel, shape = TumorType)) +
  geom_point(size=3,aes(shape=factor(TumorType))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 24, colour = "black"),
        axis.text = element_text(size = 24, colour = "black"),
        legend.text = element_text(size = 15, colour = "black"),
        legend.title = element_text(size = 20, colour = "black"))+
  geom_point(size=5) +
  scale_shape_manual(values=1:15) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  scale_color_manual(values = c(P1="black", P10="grey",
                                P2="red",P3="yellow",
                                P4="green",P5="blue4",
                                P6="pink",P7="sienna",
                                P8="magenta",P9="orange",
                                P11="royalblue1",P12="olivedrab",
                                P13="khaki",P14="palegreen2",
                                P15="orchid",P16="tan",
                                P17="cadetblue2",P18="cornsilk",
                                P19="hotpink2",P20="gold2")
  )

dev.off()


##################################volcano plot 
##organize dataframe
rownames(resOrdered)
resOrdered$SYMBOL <- rownames(resOrdered)
out <- data.frame(resOrdered)
library(tidyr)
out<-out %>% drop_na()
###volcano plot
out$Diffexpressed<-"Not Significant"
out$Diffexpressed[out$log2FoldChange > 0.5 & out$pvalue < 0.05] <- "Up"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
out$Diffexpressed[out$log2FoldChange < -0.5 & out$pvalue < 0.05] <- "Down"
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
head(out)
out$delabel <- NA
up <- rownames(out[out$Diffexpressed=="Up", ])[1:20]
down <- rownames(out[out$Diffexpressed=="Down", ])[1:20]

##label name
delabel<- c("CLDN4","CLDN7","ERBB3","CDH1","SPINT2","TFAP2A","PTPRF","MUC1","TACSTD2","F11R","KRT8","DSP","FOXA1","CD24",
"GATA3","MAL2","PERP","TFF3","LTF",
"BGN","COL1A2","IGFBP7","POSTN","COL3A1","AEBP1","IGHG3","LAPTM5","IGKC","IGHA1","MMP2","SPARC","COL5A1","APOE",
"A2M","C1QC","CD68","FN1","FBN1","VIM","CCL2","CXCL9","SPP1")

delabel<- c("CHST1","CISD3","AGR2","MLLT6","PIP4K2B","ANKRD17","LASP1","PSCA","PCGF2","TBC1D3C","JPT1",
            "PSMB3","RPL23","LLGL2", "CANX","APOD","H3C13","PSMD3","TFF3", "MUCL1","LTF",
            "KRT81", "IGFBP4","S100A1","COL6A1","KLK5","KRT5","C19orf33","CSTA","HLA-A","UGCG",
            "MFGE8","SFRP1","KRT14","SERPINA3","CCND1","FGFR1","KLK8","MAL2","CCL28", "EGR1")

delabel<- c(up,down) 
legend_title<-""
#out$delabel[out$diffexpressed != "NO"] <- out$SYMBOL[out$diffexpressed != "NO"]
library(ggrepel)
pdf("./Figures/NonmetvsPrimary_volcanoplot.pdf", width = 9, height = 6)
#pdf("./Figures/TumorvsStroma_volcanoplot.pdf", width = 9, height = 6)
ggplot(data=out, aes(x=log2FoldChange, y=-log10(pvalue), col=Diffexpressed)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel(
    data = out[tolower(out$SYMBOL)%in%tolower(delabel),], #subset(out, -log10(out$padj)>15)
    aes(label = SYMBOL), size=5, color="red4", box.padding = unit(0.35, "lines"),
    max.overlaps =40,
    fontface = "italic",
    point.padding = unit(0.4, "lines")) +
  geom_text(aes(label=ifelse(log10(padj)>15,as.character(SYMBOL),'')),hjust=0,vjust=0 ) +
  ylim(0,10) + ##for LN vs primary
  xlab("Log2FoldChange") +
  ylab("P value (-log10)") +
  scale_color_manual(legend_title, values = c("blue", "grey","red")) +
  theme_bw() +
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = 'black'),
    axis.text = element_text(size=18, color = "black"), 
    axis.title = element_text(size=18),
    legend.text = element_text(size=16),
    legend.title = element_text(size=13),
    legend.position = c(0.9,0.90))


dev.off()
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

fil_lst <- c("asthma", "allograft_rejection","autoimmune_thyroid_disease","ribosome",
             "type_i_diabetes_mellitus", "viral_myocarditis","graft_versus_host_disease", 
             "leishmania_infection","parkinsons_disease","glioma","pathways_in_cancer",
             "huntingtons_disease", "small_cell_lung_cancer","non_small_cell_lung_cancer",
             "prion_diseases","chronic_myeloid_leukemia","pancreatic_cancer",
             "endometrial_cancer","alzheimers_disease","colorectal_cancer","acute_myeloid_leukemia") ## for IM

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
  theme(axis.text.x = element_text(color="black", size=14, angle=0),
        axis.text.y = element_text(color="black", size=14, angle=0)) +
  theme( axis.line = element_line(colour = "black", size = 1, linetype = "solid"))+
  theme(axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"))+ 
  theme(legend.title = element_text(colour="black", size=14))+
  theme(legend.text = element_text(colour="black", size=14),
        legend.position = c(0.7,0.35))+
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

#pdf("./Figures/LNvsPrimary_GSEA.pdf", width=10, height=7)
pdf("./Figures/NonMetvsprimary.pdf", width=10, height=10)
p1
dev.off()




