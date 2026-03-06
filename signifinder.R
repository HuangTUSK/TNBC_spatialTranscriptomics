
library(signifinder)
mtx <- read.table("TPMnormalized.txt", sep = "\t", header = T)
head(mtx)
rownames(mtx) <- mtx$X
mtx$X <- NULL

res <- ICBResponseSign(mtx, nametype = "SYMBOL", whichAssay = "norm_expr")
ov_res <- consensusOVSign(mtx, nametype = "SYMBOL", whichAssay = "norm_expr")
ov_res@colData
df1 <- data.frame(ov_res@colData)
correlationSignPlot(ov_res)
df1$Sample <- rownames(df1)

######meta data
TNBC <- read.table("./Deseq2/samplelist.txt", sep = "\t" , header = T)
TNBC$AOI <- gsub("\\..*", "", TNBC$AOI)
TNBC$Label <- gsub("\\_.*", "", TNBC$Label)


TNBC$merge <- paste(TNBC$TNBCtype, TNBC$X, sep="_")
all <- read.table("./TNBCmetaData.txt", sep = "\t", header = T) 
head(all)
all$AOISurfaceArea <- gsub("\\..*", "", all$AOISurfaceArea)
sel <- all[,c("SegmentDisplayName","AOISurfaceArea")]

mergedDF <- merge(sel, TNBC, by.x = "AOISurfaceArea", by.y = "AOI")
mergedDF$TumorType <- gsub("-","_", mergedDF$TumorType)
mergedDF$TumorType <- gsub("\\+","_P", mergedDF$TumorType)
rownames(mergedDF) <- mergedDF$SegmentDisplayName

metadf <- mergedDF[,c("MetType","Label","P53","X","type")]
#rownames(metadf) <- gsub("\\|","", rownames(metadf))
#rownames(metadf) <- gsub(" ","", rownames(metadf))
rownames(metadf) <- metadf$segment
metadf2 <- metadf[,c("Max_Column","TNBCtype","Label","SpatialType", "Tissue")]

rownames(metadf) <- gsub("\\|","", rownames(metadf))
rownames(metadf) <- gsub(" ","", rownames(metadf))

#df <- df[order(df$X),]

metadf <- metadf[metadf$type=="Tumor",]
tdf2 <- tdf[,rownames(metadf)]

CIN_res <- CINSign(dataset = mtx, inputType = "rnaseq",nametype = "SYMBOL", whichAssay = "norm_expr")
CIN_res@colData

df2 <- data.frame(CIN_res@colData)
df2$Sample <- rownames(df2)
##
cellcycle_L <- cellCycleSign(dataset = mtx,author = "Lundberg", inputType = "rnaseq",nametype = "SYMBOL", whichAssay = "norm_expr")
cellcycle_L@colData

df3 <- data.frame(cellcycle_L@colData)
df3$Sample <- rownames(df3)

cellcycle_D <- cellCycleSign(dataset = mtx,author = "Davoli", inputType = "rnaseq",nametype = "SYMBOL", whichAssay = "norm_expr")
cellcycle_D@colData

df4 <- data.frame(cellcycle_D@colData)
df4$Sample <- rownames(df4)

availableSignatures(topic = "cell cycle")

##
mitoI_res <- mitoticIndexSign(dataset = mtx, nametype = "SYMBOL", whichAssay = "norm_expr")
mitoI_res@colData

df5 <- data.frame(mitoI_res@colData)
df5$Sample <- rownames(df5)

##
DNAre_res <- DNArepSign(dataset = mtx, inputType = "rnaseq",nametype = "SYMBOL", whichAssay = "norm_expr")
DNAre_res@colData

df6 <- data.frame(DNAre_res@colData)
df6$Sample <- rownames(df6)

##
IPSOV_res <- IPSOVSign(dataset = mtx, inputType = "rnaseq",nametype = "SYMBOL", whichAssay = "norm_expr")
IPSOV_res@colData

df7 <- data.frame(IPSOV_res@colData)
df7$Sample <- rownames(df7)

##
IPS_res <- IPSSign(dataset = mtx, nametype = "SYMBOL", whichAssay = "norm_expr",hgReference = "hg38")
IPS_res@colData

df8 <- data.frame(IPS_res@colData)
df8$Sample <- rownames(df8)

##
che_res <- chemokineSign(dataset = mtx, inputType = "rnaseq",nametype = "SYMBOL", whichAssay = "norm_expr")
che_res@colData

df9 <- data.frame(che_res@colData)
df9$Sample <- rownames(df9)

##
ECM_res <- ECMSign(dataset = mtx, nametype = "SYMBOL", whichAssay = "norm_expr")
ECM_res@colData

df10 <- data.frame(ECM_res@colData)
df10$Sample <- rownames(df10)

##
matri_res <- matrisomeSign(dataset = mtx, nametype = "SYMBOL", whichAssay = "norm_expr")
matri_res@colData

df11 <- data.frame(matri_res@colData)
df11$Sample <- rownames(df11)

##
lip_res <- lipidMetabolismSign(dataset = mtx, nametype = "SYMBOL", whichAssay = "norm_expr")
lip_res@colData

df12 <- data.frame(lip_res@colData)
df12$Sample <- rownames(df12)

##
availableSignatures(topic = "epithelial")
ETM_res <- EMTSign(dataset = mtx, author = "Mak",nametype = "SYMBOL", whichAssay = "norm_expr",hgReference = "hg38")
ETM_res@colData

df13 <- data.frame(ETM_res@colData)
df13$Sample <- rownames(df13)

##
ExImm_res <- expandedImmuneSign(dataset = mtx, nametype = "SYMBOL", whichAssay = "norm_expr")
ExImm_res@colData

df14 <- data.frame(ExImm_res@colData)
df14$Sample <- rownames(df14)

##
IFN_res <- IFNSign(dataset = mtx, nametype = "SYMBOL", whichAssay = "norm_expr")
IFN_res@colData

df15 <- data.frame(IFN_res@colData)
df15$Sample <- rownames(df15)

##
availableSignatures(topic = "IPS")
Iscore_H <- immunoScoreSign(dataset = mtx, author = "Hao",inputType = "rnaseq", nametype = "SYMBOL", whichAssay = "norm_expr",hgReference = "hg38")
Iscore_H@colData

df16 <- data.frame(Iscore_H@colData)
df16$Sample <- rownames(df16)

##
Iscore_R <- immunoScoreSign(dataset = mtx, author = "Roh", inputType = "rnaseq", nametype = "SYMBOL", whichAssay = "norm_expr",hgReference = "hg38")
Iscore_R@colData

df17 <- data.frame(Iscore_R@colData)
df17$Sample <- rownames(df17)

##
Icyt_H <- immuneCytSign(dataset = mtx, author = "Rooney", inputType = "rnaseq", nametype = "SYMBOL", whichAssay = "norm_expr",hgReference = "hg38")
Icyt_H@colData

df18 <- data.frame(Icyt_H@colData)
df18$Sample <- rownames(df18)

##
Icyt_D <- immuneCytSign(dataset = mtx, author = "Davoli",inputType = "rnaseq", nametype = "SYMBOL", whichAssay = "norm_expr",hgReference = "hg38")
Icyt_D@colData

df19 <- data.frame(Icyt_D@colData)
df19$Sample <- rownames(df19)

##
Tinflam <- TinflamSign(dataset = mtx, author = "Ayers", nametype = "SYMBOL", whichAssay = "norm_expr",hgReference = "hg38")
Tinflam@colData

df20 <- data.frame(Tinflam@colData)
df20$Sample <- rownames(df20)

##
TGFB_res <- TGFBSign(dataset = mtx, nametype = "SYMBOL", whichAssay = "norm_expr")
TGFB_res@colData

df21 <- data.frame(TGFB_res@colData)
df21$Sample <- rownames(df21)

##
VEGF_res <- VEGFSign(dataset = mtx, nametype = "SYMBOL", whichAssay = "norm_expr")
VEGF_res@colData

df22 <- data.frame(VEGF_res@colData)
df22$Sample <- rownames(df22)

##
ADO_res <- ADOSign(dataset = mtx, nametype = "SYMBOL", whichAssay = "norm_expr")
ADO_res@colData

df23 <- data.frame(ADO_res@colData)
df23$Sample <- rownames(df23)

##
APM_res <- APMSign(dataset = mtx, author = "Wang",inputType = "rnaseq", nametype = "SYMBOL", whichAssay = "norm_expr",hgReference = "hg38")
APM_res@colData

df24 <- data.frame(APM_res@colData)
df24$Sample <- rownames(df24)

##
Auto_res <- autophagySign(dataset = mtx, author = "Xu", nametype = "SYMBOL", whichAssay = "norm_expr",hgReference = "hg38")
Auto_res@colData

df25 <- data.frame(Auto_res@colData)
df25$Sample <- rownames(df25)

###
com_res <- CombinedSign(dataset = mtx,nametype = "SYMBOL", whichAssay = "norm_expr",hgReference = "hg38")
com_res@colData

df26 <- data.frame(com_res@colData)
df26$Sample <- rownames(df26)

###
gly_res <- glycolysisSign(dataset = mtx,nametype = "SYMBOL",author = "Zhang", whichAssay = "norm_expr")
gly_res@colData

df27 <- data.frame(gly_res@colData)
df27$Sample <- rownames(df27)

##
HRD_res <- HRDSSign(dataset = mtx,nametype = "SYMBOL", whichAssay = "norm_expr")
HRD_res@colData

df29 <- data.frame(HRD_res@colData)
df29$Sample <- rownames(df29)

##
Hypoxia_res <- hypoxiaSign(dataset = mtx,nametype = "SYMBOL", whichAssay = "norm_expr")
Hypoxia_res@colData

df30 <- data.frame(Hypoxia_res@colData)
df30$Sample <- rownames(df30)

##
ICBR_res <- ICBResponseSign(dataset = mtx,nametype = "SYMBOL", whichAssay = "norm_expr")
ICBR_res@colData

df31 <- data.frame(ICBR_res@colData)
df31$Sample <- rownames(df31)

##

###
lst <- list (df2, df3, df4, df5,df6, df7,df8, df9, df10,
            df11, df12, df13, df14, df15,
            df16, df17, df18, df19, df20, df21, df22,
            df23,df24,df25,df26,df27,df29,df30,
            df31)
matx <- purrr::reduce(.x = lst, merge, by = 'Sample', all = T)
matx$IPS_Charoentong <- NULL

write.table(matx, "Signifinder_matrix.txt", sep = "\t", quote = F, row.names = F)


#####################correlation 
matx <- read.table("Signifinder_matrix.txt", sep = "\t", header = T)
library(reshape2)
matx$segment <- gsub("\\.", "", matx$Sample)

#metadf <- allsub_merged
####meta data 
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
mergedDF$SegmentDisplayName <- gsub("\\|","", mergedDF$SegmentDisplayName)
mergedDF$SegmentDisplayName <- gsub(" ","", mergedDF$SegmentDisplayName)
mergedDF$SegmentDisplayName 
merge_sigMatx <- merge(mergedDF,matx,by.x="SegmentDisplayName",by.y="segment")
merge_sigMatx$SpatialType <- merge_sigMatx$X

colnames(merge_sigMatx)
rownames(merge_sigMatx) <- merge_sigMatx$SegmentDisplayName

###select subgroup
merge_sigMatx_sub <- merge_sigMatx[merge_sigMatx$TNBCtype=="UNS",]

sigM_cor <- merge_sigMatx_sub[,c(23:48)]
cor_matrix <- cor(sigM_cor)
# Reshape the correlation matrix for ggplot
cor_melted <- melt(cor_matrix)
ggplot(cor_melted, aes(Var1, Var2, fill = value)) + 
  geom_tile() + 
  scale_fill_gradient2(midpoint = 0, low = "red", high = "blue", mid = "white") + 
  theme_minimal() + 
  labs(title = "Correlation Heatmap", x = "Variables", y = "Variables") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

###heatmap 
label_matrix <- cor_matrix*100
pheatmap(
  cor_matrix,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "ward.D2",
 # annotation_col = data.frame(Group = factor(groups)),
 # annotation_row = data.frame(Group = factor(groups)),
  color = colorRampPalette(c("blue", "white", "red"))(100),
  fontsize_row = 6,
  fontsize_col = 6,
  legend = TRUE
 #display_numbers = label_matrix, fontsize_number=6,
)


#####heatmap
nonlst <- c("ImmunoScore_Hao", "ECM_Chakravarthy_up","ECM_Chakravarthy_down",
            "CellCycle_Lundberg","IPSOV_Shen","ICBResponse_Chen_responder",
            "ICBResponse_Chen_nonresponder","Combined_Thompson","ImmuneCyt_Rooney")
matx <- matx[,!names(matx) %in% nonlst]
matx$CellCycle_Lundberg
rownames(matx) <- matx$Sample
matx$Sample <- NULL

Lmatx <- apply(matx, 1, function(x) (x - mean(x)) / sd(x))
#Lmatx <- apply(matx, 1, function(x) x / max(x))
#tdf <- data.frame(t(Lmatx))
#tdf$ECM_Chakravarthy_down
pheatmap(Lmatx,
         col=colors,
         show_colnames=F,
         cluster_rows=T)

#####annotation 
metadf <- allsub_merged
#metadf <- df ## df from spatial Deseq2 

#rownames(metadf) <- gsub("\\|","", rownames(metadf))
#rownames(metadf) <- gsub(" ","", rownames(metadf))
rownames(metadf) <- metadf$segment
#metadf2 <- metadf[,c("TNBCtype","Label","SpatialType")]
#df <- df[order(df$X),]

####select DIF & Exclused
metadf <- metadf[metadf$type=="Tumor" & metadf$Tissue=="LN",] #breast
#metadf <- metadf[metadf$Max_Column=="ConsensusOV_Chen_DIF",] ##& metadf$SpatialType=="Excluded"
colnames(Lmatx) <- gsub("\\.", "", colnames(Lmatx))
Lmatx2 <- Lmatx[,metadf$segment]

metadf$TNBCtype > metadf$TNBCtype.x
metadf$SpatialType <- metadf$X
colnames(metadf)
metadf2 <- metadf[,c("TNBCtype.x","Label","SpatialType","Tissue")]

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
  SpatialType = c(Excluded = "olivedrab",
                 Ignored = "tan",
                 Inflamed = "royalblue"),
  Tissue = c(breast = "black",
             LN = "grey"),
  TNBCtype = c(BL1="gold",
               BL2="cadetblue",
               IM="pink",
               LAR="palegreen",
               M="blue4",
               MSL="grey51",
               UNS="orange4"))

metadf2 <- metadf2[order(metadf$TNBCtype),]
Lmatx2 <- Lmatx2[,rownames(metadf2)]
pheatmap(Lmatx2, cluster_rows=T, show_rownames=T,
         show_colnames=F, cluster_cols=T, fontsize_row = 10,
         annotation_col=metadf2,
         annotation_colors = ann_col)



