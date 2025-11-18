
cyber <- read.table("CIBERSORTx_Job82_Adjusted.txt", sep = "\t", header = T)
head(cyber)
rowSums(cyber[,2:23])

TNBC <- read.table("samplelist.txt", sep = "\t" , header = T)
head(TNBC)
TNBC$DSPID
TNBC$AOI <- gsub("\\..*", "", TNBC$AOI)

all <- read.table("./TNBCmetaData.txt", sep = "\t", header = T) 
head(all)
all$AOISurfaceArea <- gsub("\\..*", "", all$AOISurfaceArea)

sel <- all[,c("SegmentDisplayName","AOISurfaceArea")]

mergedDF <- merge(sel, TNBC, by.x = "AOISurfaceArea", by.y = "AOI")

mergecyber <- merge(mergedDF, cyber, by.x = "SegmentDisplayName", by.y = "Mixture")
nrow(mergecyber)
head(mergecyber)

tumor <- mergecyber[mergecyber$type=="Tumor",]
stroma <- mergecyber[mergecyber$type=="Stroma",]

colnames(tumor)
resDF<-data.frame()
for (column in colnames(mergecyber)[c(24:45)]) {
  resDF<-rbind (
    resDF,
    data.frame(Signature=column,
               pvalue=wilcox.test(mergecyber[mergecyber$type=="Tumor",][,column], mergecyber[mergecyber$type=="Stroma",][,column])$p.value,
               Tumormean= mean(mergecyber[mergecyber$type=="Tumor",][,column], na.rm=TRUE),
               Stromamean= mean(mergecyber[mergecyber$type=="Stroma",][,column], na.rm=TRUE)
    )
  )
}

resDF

#####location
tumor <- mergecyber[mergecyber$type=="Tumor",]
stroma <- mergecyber[mergecyber$type=="Stroma",]

tumor$MetType
tumor <- tumor[1:39,]
##tumor
resDF<-data.frame()
for (column in colnames(tumor)[c(24:45)]) {
  resDF<-rbind (
    resDF,
    data.frame(Signature=column,
               pvalue=wilcox.test(tumor[tumor$Tissue=="breast",][,column], 
                                  tumor[tumor$Tissue=="LN",][,column])$p.value,
               Primarymean= mean(tumor[tumor$Tissue=="breast",][,column], na.rm=TRUE),
               LNmean= mean(tumor[tumor$Tissue=="LN",][,column], na.rm=TRUE)
    )
  )
}

resDF
write.table(resDF, "Primary_LN_Immunecell.txt", quote = F, sep = "\t", row.names = F)

######
df <- na.omit(resDF)
df$EnrichedScore <- ifelse(df$Primarymean > df$LNmean,-(df$Primarymean), df$LNmean)
df$Significant <- ifelse(df$pvalue < 0.05, "Yes", "No")
df$Annotation_Color = c("#4CAF50", "#C5E1A5", "#8BC34A", "#1B5E20", "#795548",
                        "tan2","cadetblue1","#3c5488","tan","cadetblue4","royalblue",
                        "#e64b35", "#4dbbd5","#00a087",  "#f39b7f","wheat4","hotpink3",
                        "slategray4","sienna","magenta","khaki","palegreen2") # Custom colors for annotation


p1 <- ggplot(df, aes(x = reorder(Signature, EnrichedScore), y = EnrichedScore, fill = Significant)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Flip axes
  scale_fill_manual(values = c("Yes" = "black", "No" = "gray")) + 
  labs(x = NULL, y = "Primary enriched | LN enriched", fill = "Significant\n(p < 0.05)") +
  theme_minimal() +
  theme(#axis.text.y = element_blank(), 
    axis.text.y = element_text(size = 10, color = "black"),
    axis.ticks.y = element_blank())

p1

p2 <- ggplot(df, aes(x = reorder(Signature, EnrichedScore), y = 1, fill = Annotation_Color)) +
  geom_tile() +
  coord_flip() +
  scale_fill_identity() +  # Use custom colors directly
  # theme_void() 
  theme(
    axis.title.y = element_blank(),   # Remove y-axis title
    axis.text.y = element_text(size = 10, color = "black"),  # Show y-axis text
    #   axis.ticks.y = element_blank(),   # Hide y-axis ticks if not needed
    #  axis.line.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x =element_blank()# Hide y-axis line if not needed
  )
p2
final_plot <- ggarrange(p2, p1, ncol = 2, widths = c(0.3, 1))

pdf("./Figures/Cybersort_LN_nonMetPrimarybarplot.pdf", width = 6, height = 5)
p1
dev.off()

###stroma 

resDF<-data.frame()
for (column in colnames(stroma)[c(25:45)]) {
  resDF<-rbind (
    resDF,
    data.frame(Signature=column,
               pvalue=wilcox.test(stroma[stroma$Tissue=="breast",][,column], 
                                  stroma[stroma$Tissue=="LN",][,column])$p.value,
               Breastmean= mean(stroma[stroma$Tissue=="breast",][,column], na.rm=TRUE),
               Lynmean= mean(tumor[stroma$Tissue=="LN",][,column], na.rm=TRUE)
    )
  )
}

resDF

####plot 
head(tumor)
colnames(tumor)
tumor[,c("SegmentDisplayName","Patient")]
stroma[,c("SegmentDisplayName","Patient")]
library(tidyr)
library(tidyverse)

colnames(mergecyber)
mergecyber$MetType
mergecyber2 <- mergecyber

for (i in 24:45) {
  mergecyber2[[i]] <- mergecyber2[[i]] / mergecyber2[[49]]
}

mergecyber2$MetType
#mergecyber2 <- mergecyber2[1:78,]. ##for no nonmet
mergecyber2$Absolute.score..sig.score.
df_new <- tidyr::gather(mergecyber2,key='celltype',value='Fraction',24:45)
df_new2 <- df_new[df_new$type=="Tumor",]
df_new$patient_type <- paste(df_new$Patient, df_new$type, sep = "")
df_new2$MetType <- sub("\\..*","", df_new2$MetType)
df_new2$MetType <- gsub('[0-9]+', '', df_new2$MetType)
df_new2$MetType <- gsub(' ', '', df_new2$MetType)
##plot
pdf("./Figures/Cybersort_stackedplot_mettype.pdf", width = 16, height = 5.5)
ggplot(df_new2, aes(fill=celltype, y=Fraction, x=DSPID)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("Primary VS. Ln tumor") +
  xlab("") +
  ylab("Cell components (Relative value)") +
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
  facet_grid(.~MetType, scales = "free_x") +
  scale_fill_manual(values = c("hotpink","hotpink3","khaki1","khaki3",
                               "skyblue","skyblue2","skyblue4",
                               "sienna1","sienna3",
                               "#00a087",
                               "#3c5488", 
                               "#f39b7f",
                               "darkseagreen1","darkseagreen3",
                               "thistle",
                               "azure1","azure2","azure3","azure4",
                               "slategray1","slategray2","slategray3"
                               ))

dev.off()

#########all
#df_new <- tidyr::gather(mergecyber2,key='celltype',value='Fraction',24:45)
#df_new2 <- df_new[df_new$type=="Tumor",]
#df_new$patient_type <- paste(df_new$Patient, df_new$type, sep = "")

pdf("./Figures/Cybersort_all_stackedplot_mettype.pdf", width = 16, height = 5.5)
ggplot(df_new2, aes(fill=celltype, y=Fraction, x=DSPID)) + 
  geom_bar(position="stack", stat="identity") +
 # ggtitle("Primary VS. Ln tumor") +
  xlab("") +
  ylab("Cell components (Relative value)") +
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
  facet_grid(.~MetType, scales = "free_x") +
  scale_fill_manual(values = c("hotpink","hotpink3","khaki1","khaki3",
                               "skyblue","skyblue2","skyblue4",
                               "sienna1","sienna3",
                               "#00a087",
                               "#3c5488", 
                               "#f39b7f",
                               "darkseagreen1","darkseagreen3",
                               "thistle",
                               "azure1","azure2","azure3","azure4",
                               "slategray1","slategray2","slategray3"
  ))

dev.off()



###tumor 
tumor_set <- df_new[df_new$type=="Tumor",]
ggplot(tumor_set, aes(fill=celltype, y=Fraction, x=DSPID)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("Tumor") +
  xlab("") +
  ylab("Immune Cell Components") +
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
  facet_grid(.~patient_type, scales = "free_x") +
  scale_fill_manual(values = c("hotpink","hotpink3","khaki1","khaki3",
                               "skyblue","skyblue2","skyblue4",
                               "sienna1","sienna3",
                               "#00a087",
                               "#3c5488", 
                               "#f39b7f",
                               "darkseagreen1","darkseagreen3",
                               "thistle",
                               "azure1","azure2","azure3","azure4",
                               "slategray1","slategray2","slategray3"
  ))

##stroma
tumor_set <- df_new[df_new$type=="Stroma",]
ggplot(tumor_set, aes(fill=celltype, y=Fraction, x=DSPID)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("Stroma") +
  xlab("") +
  ylab("Immune Cell Components") +
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
  facet_grid(.~patient_type, scales = "free_x") +
  scale_fill_manual(values = c("hotpink","hotpink3","khaki1","khaki3",
                               "skyblue","skyblue2","skyblue4",
                               "sienna1","sienna3",
                               "#00a087",
                               "#3c5488", 
                               "#f39b7f",
                               "darkseagreen1","darkseagreen3",
                               "thistle",
                               "azure1","azure2","azure3","azure4",
                               "slategray1","slategray2","slategray3"
  ))


##############compare among subtypes 

head(mergecyber)
colnames(mergecyber)
head(subtypeMerge)
subtypeMerge$Genefu <- subtypeMerge$subtype.x
subtypeMerge$TNBCtype <- subtypeMerge$subtype.y
colnames(subtypeMerge)
subtype <- subtypeMerge[,c(1,19,20,21,22,23,31,32)]


allmerge <- merge(subtype, mergecyber, by.x= "SegmentDisplayName", by.y = "SegmentDisplayName")

#write.table(allmerge, "Tumor_subtype_immunecellFraction.txt", sep = "\t", quote = F)
allmerge <- read.table("Tumor_subtype_immunecellFraction.txt",sep = "\t",header = T)
####
colnames(mergecyber)
df_new <- tidyr::gather(allmerge,key='cellType',value='Fraction',25:46,-1)

resDF<-data.frame()
for (column in colnames(mergecyber)[c(24:45)]) {
  resDF<-rbind (
    resDF,
    data.frame(Signature=column,
               summary(aov(mergecyber[,column] ~ mergecyber[,21]))[[1]][["Pr(>F)"]][1],
              # pvalue=wilcox.test(allmerge[allmerge$Tissue=="breast",][,column], 
              #                    allmerge[allmerge$Tissue=="LN",][,column])$p.value,
              # BL1tmean= mean(mergecyber[mergecyber$TNBCtype=="BL1",][,column], na.rm=TRUE),
              # BL2nmean= mean(mergecyber[mergecyber$TNBCtype=="BL2",][,column], na.rm=TRUE),
              # IMnmean= mean(mergecyber[mergecyber$TNBCtype=="IM",][,column], na.rm=TRUE),
              # LARnmean= mean(mergecyber[mergecyber$TNBCtype=="LAR",][,column], na.rm=TRUE),
              # MSLnmean= mean(mergecyber[mergecyber$TNBCtype=="MSL",][,column], na.rm=TRUE),
              # UNSnmean= mean(mergecyber[mergecyber$TNBCtype=="UNS",][,column], na.rm=TRUE)
              ExcludedMean=mean(mergecyber[mergecyber$X=="Excluded",][,column], na.rm=TRUE),
              IgnoredMean=mean(mergecyber[mergecyber$X=="Ignored",][,column], na.rm=TRUE),
              InflamedMean=mean(mergecyber[mergecyber$X=="Inflamed",][,column], na.rm=TRUE)
              
    )
  )
  
  
}

resDF
#write.table(resDF, "TNBCtype_Immunecell.txt", sep = "\t", quote = F, row.names = F)
write.table(resDF, "Spatialtype_Immunecell.txt", sep = "\t", quote = F, row.names = F)

##plot
library(ggsignif)
library(ggplot2)
library(tidyverse)
library(ggpubr)


p1<-ggplot(data = df_new, aes(x=cellType, y=Fraction)) +
  geom_boxplot(aes(fill=Genefu), outlier.size = 0.1) +
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
p1 + stat_compare_means(aes(group = Genefu),method="anova", label = "p.signif",
                        label.y = 0.38, size =6, hide.ns = TRUE) + ggtitle("")



###
colnames(mergecyber)
df_new <- tidyr::gather(mergecyber,key='cellType',value='Fraction',18:39,-1)

p1<-ggplot(data = df_new, aes(x=cellType, y=Fraction)) +
  geom_boxplot(aes(fill=type), outlier.size = 0.1) +
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
p1 + stat_compare_means(aes(group = type),method="wilcox.test", label = "p.signif",
                        label.y = 0.6, size =6, hide.ns = TRUE) + ggtitle("")


####################subtype 
colnames(allmerge)


p1<-ggplot(data = allmerge, aes(x=TNBCtype, y=TumorSize)) +
  geom_boxplot(aes(fill=TNBCtype), outlier.size = 0.1) +
  ylim(0,8) +
  labs(fill = "") +
  # ylab("Expression (log10)") +
  ylab("Tumor size")+
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
        axis.text.x = element_text(size=14, color = "black", angle = 45, hjust=1),
        # axis.text.x = element_text(size=15, color = "black"), axis.title = element_text(size=22),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        plot.title = element_text(size = 27, hjust = 0.5),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        #legend.title = element_text(size=13),
        #legend.position = c(0.9,0.7)
  ) #c(0.9,0.8)
p1 
p1 + stat_compare_means(aes(group = TNBCtype),method="anova", label = "p.signif",
                        label.x = 0.5, label.y = 7.5, size =8, hide.ns = TRUE) + ggtitle("")






