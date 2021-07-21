#### Set Up Libraries ####

#biocLite("Gviz")
#biocLite("org.Hs.eg.db")
library("here")
library("DESeq2")
library("dplyr")
library("tidyr")
library("AnnotationDbi")
library("KEGG.db")
library("RColorBrewer")
library("pheatmap")
library("ggplot2")
library("ggrepel")
library("extrafont")
library("scales")
library("vsn")
library("biomaRt")
library("categoryCompare")
library("RCy3")
library("GOstats")
library("tximport")
library("readr")
library("ReportingTools")
#BiocManager::install("ReportingTools")


#### Load files ####
sampleTable<-data.frame(SampleID=c("SE38-16_S6",
                                   "SE38-19_S8","SE38-20_S9","SE38-25_S10",
                                   "SE38-29_S11","SE38-30_S12"),
                        Genotype=c("WT","KO","KO","KO","WT","WT"))

sampleTable<-sampleTable[sampleTable$SampleID!="SE38-13_S5"&sampleTable$SampleID!="SE38-17_S7",] 
#samples 13 and 17 are not the same as the others

#### Establish file paths ####
dir<-"./quantsSE"
samps<-list.dirs(dir,full.names=F,recursive=F)
files <- file.path(dir, samps, "quant.sf")
names(files) <- gsub("_quant","",samps)
setdiff(sampleTable$SampleID,names(files))
files<-files[names(files)%in%sampleTable$SampleID]
all(file.exists(files)) 
files[!file.exists(files)]
row.names(sampleTable)<-sampleTable$SampleID
idx<-match(rownames(sampleTable),names(files))
files<-files[idx]
identical(names(files),rownames(sampleTable))
sampleTable<-sampleTable[sampleTable$SampleID%in%names(files),]

quant<-read.table(files[1],sep="\t",header=T,as.is=T,check.names = F)
txs<-quant$Name
txs<-sapply(strsplit(txs,"[.]"),"[",1)


ensembl<-useEnsembl("ensembl",dataset="mmusculus_gene_ensembl")
tx2gene<-getBM(attributes=c('ensembl_transcript_id',"ensembl_gene_id"),
               filters="ensembl_transcript_id",values=txs,mart=ensembl)
#rm(quant,txs)

#### Run tximport ####
txi <- tximport(files, type="salmon", tx2gene=tx2gene, ignoreTxVersion = T, #countsFromAbundance=c("lengthScaledTPM"),
                dropInfReps = T)
save(txi,file="txi.RData")
ddsTxi <- DESeqDataSetFromTximport(txi,colData = sampleTable,design = ~ Genotype)
ddsTxi <- estimateSizeFactors(ddsTxi,type="poscounts")
ddsTxi <- DESeq(ddsTxi)
save(txi,ddsTxi,file="txi.RData")

#### Create DESeq objects ####
#filtering
ddsTxi<- ddsTxi[which(mcols(ddsTxi)$betaConv),]
keep<-apply(counts(ddsTxi),1,function(x) length(x[x>5])>=2) 
sum(keep) #sum(keep) is 17324
ddsTxi <- ddsTxi[keep,]
#ddsTxi <- ddsTxi[!row.names(counts(ddsTxi))%in%acinarGenes,]

#dds<-DESeq(dds,modelMatrixType = "expanded",betaPrior = T)
resultsNames(ddsTxi)

#### Get annotation for genes from Ensembl ####
#(uses biomaRt library)
# Open connection to Ensembl           hsapiens_gene_ensembl
#ensembl<-useEnsembl("ensembl",dataset="mmusculus_gene_ensembl")

# Get annotation (Ensembl ID, Gene Symbol, gene biotype) for each ID in the featureCounts table
#filteredIDs<-sapply(strsplit(rownames(fc_filtered),"[.]"),"[",1)
filteredIDs<-rownames(counts(ddsTxi))
genemap<-getBM(attributes=c("ensembl_gene_id","external_gene_name","gene_biotype","entrezgene"),
               filters="ensembl_gene_id",values=filteredIDs,mart=ensembl)

# Make index linking rows from featureCounts to rows in genemap
genemap_idx<-match(filteredIDs,genemap$ensembl_gene_id)


counts<-counts(ddsTxi,normalized=T)
symbols<-data.frame("symbols"=genemap$external_gene_name[genemap_idx])
row.names(symbols)<-row.names(counts)
counts<-data.frame(genemap$external_gene_name[genemap_idx],counts,check.names=F)
names(counts)[1]<-"Symbols"
head(counts)
write.csv(counts,file="counts.csv",row.names=T,quote=F)

#### Pick normalization method ####

ntd <- normTransform(ddsTxi)
rld <- rlog(ddsTxi, blind=FALSE) #slow, too slow for large datasets
# blind = F because differences in groups should not be interpreted as noise
# we expect to see many differences in genes due to experimental groups
vsd <- varianceStabilizingTransformation(ddsTxi, blind=FALSE)
vsd.fast <- vst(ddsTxi, blind=FALSE)

df <- bind_rows(
  as_tibble(log2(counts(ddsTxi, normalized=TRUE)+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(rld)) %>% mutate(transformation = "rlog"),
  as_tibble(assay(vsd)) %>% mutate(transformation = "vst"),
  as_tibble(assay(ntd)) %>% mutate(transformation = "ntd"),
  as_tibble(assay(vsd.fast)) %>% mutate(transformation = "vst.fast")
)
colnames(df)[1:2] <- c("x", "y")  
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)

#### Diagnostic MA plots to help choose transform ####
#uses vsn package
msd<-meanSdPlot(assay(ddsTxi),plot=F)
msd$gg + ggtitle("Untransformed Data") + 
  scale_fill_gradient(low = "yellow", high = "darkred") + 
  scale_y_continuous(limits = c(0, 100))
ggsave("UntransformedMSD.png",path="./Figures/QC",width = 15, height = 15, units = "cm")
msd<-meanSdPlot(assay(ntd),plot=F)
msd$gg + ggtitle("Normal Transformation") + 
  scale_fill_gradient(low = "yellow", high = "darkred") + 
  scale_y_continuous(limits = c(0, 3))
ggsave("NormTransformedMSD.png",path="./Figures/QC",width = 15, height = 15, units = "cm")
msd<-meanSdPlot(assay(rld),plot=F)
msd$gg + ggtitle("Regularized Log Transformation") + 
  scale_fill_gradient(low = "yellow", high = "darkred") + 
  scale_y_continuous(limits = c(0, 3))
ggsave("RLTransformedMSD.png",path="./Figures/QC",width = 15, height = 15, units = "cm")
msd<-meanSdPlot(assay(vsd),plot=F)
msd$gg + ggtitle("Variance Stabilizing Transformation") + 
  scale_fill_gradient(low = "yellow", high = "darkred") + 
  scale_y_continuous(limits = c(0, 3))
ggsave("VStransformedMSD.png",path="./Figures/QC",width = 15, height = 15, units = "cm")
msd<-meanSdPlot(assay(vsd.fast),plot=F)
msd$gg + ggtitle("Fast Variance Stabilizing Transformation") + 
  scale_fill_gradient(low = "yellow", high = "darkred") + 
  scale_y_continuous(limits = c(0, 3)) #this one has the flattest line
ggsave("fastVStransformedMSD.png",path="./Figures/QC",width = 15, height = 15, units = "cm")

#### Look at top 50 most expressed genes ####
#sanity check
select <- order(rowMeans(counts(ddsTxi,normalized=TRUE)),
                decreasing=TRUE)[1:50]
df <- as.data.frame(colData(ddsTxi)[,"Genotype"][order(colData(ddsTxi)$Genotype)])
names(df)<-"Genotype"
rownames(df) <- colnames(ddsTxi)[order(colData(ddsTxi)$Genotype)]
rownames(vsd.fast)<-genemap$external_gene_name[genemap_idx]

ann_colors = list(
  Genotype=c("#ff8000",
             "#000080"))

names(ann_colors[[1]])<-levels(df$Genotype)
png(filename = "./Figures/QC/Top50Heatmap.png",width=15,height=20,units="cm",res=300)
pheatmap(assay(vsd.fast)[select,order(colData(vsd.fast)$Genotype)],
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "PiYG")))(100),
         cluster_rows=FALSE, show_rownames=T,
         cluster_cols=FALSE, annotation_col=df,
         annotation_colors = ann_colors[1])
dev.off()
#rownames(vsd.fast)<-rownames(fc_filtered)

#### Diagnostic plots ####
## Make PCA plot of VST-transformed data
# All samples
pcaData <- plotPCA(vsd.fast,intgroup=c("Genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Genotype,shape=Genotype)) +
  geom_point(size=3,alpha = 0.6) +
  scale_color_manual(name="Genotype",values = ann_colors$Genotype)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: \n",percentVar[2],"% variance")) + 
  coord_fixed() +
  guides(colour=guide_legend(ncol=2))+
  theme(plot.title = element_text(family = "Arial", color="black",  size=32, hjust=0)) +
  theme(axis.title = element_text(family = "Arial", color="black", size=18))+
  theme(axis.text.x = element_text(family = "Arial",colour="black",size=14))+
  theme(axis.text.y = element_text(family = "Arial",colour="black",size=14))+
  theme(legend.text = element_text(family = "Arial",colour="black",size=14))+
  theme(legend.title = element_text(family = "Arial",colour="black",size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("PCAall.png",path="./Figures/QC",width = 20, height = 10, units = "cm")

# Plot clustered heatmap of samples
colours<-colorRampPalette(rev(brewer.pal(9,'PiYG')))(255)
sampleDists<-dist(t(assay(vsd.fast)[,order(colData(vsd.fast)$Genotype)]))
sdm<-as.matrix(sampleDists)
png(filename = "./Figures/QC/SamplesHeatmap.png",width=15,height=10,units="cm",res=300)
pheatmap(sdm, cluster_rows=FALSE, show_rownames=T,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="PiYG")))(100),
         cluster_cols=FALSE, annotation_col=df,
         annotation_colors = ann_colors)
dev.off()

#### Get Results ####
resWTvsKO <- results(ddsTxi, contrast=c("Genotype","KO","WT"))
summary(resWTvsKO)
write.csv(as.data.frame(resWTvsKO),file="results_all.csv",row.names=T,quote=F)

#### Results ####
### Get list of DE genes between Genotypes
resWTvsKO$ensembl_gene_id<-filteredIDs
resWTvsKO$external_gene_name<-genemap$external_gene_name[genemap_idx]
resWTvsKO$gene_biotype<-genemap$gene_biotype[genemap_idx]
resWTvsKO$entrezgene<-genemap$entrezgene[genemap_idx]
ind1<-which.min(resWTvsKO$padj)
g1<-resWTvsKO$external_gene_name[ind1]
ind2<-which.max(resWTvsKO$log2FoldChange)
g2<-resWTvsKO$external_gene_name[ind2]
d<-plotCounts(ddsTxi, gene=ind1, intgroup="Genotype",transform=F,
              main =g1,returnData = T)
p<-ggplot(d,aes(x=Genotype, y=count))
p+geom_boxplot(aes(fill=Genotype),alpha = 0.6)+ #geom_violin or geom_boxplot+
  geom_point(aes(fill=Genotype),colour="black",shape=21,
             position=position_jitter(w=0.05,h=0),
             na.rm=T,size=1.5)+
  guides(colour=guide_legend(ncol=2),fill=guide_legend(ncol=2))+
  scale_color_manual(name="Genotype",values = ann_colors$Genotype)+
  scale_fill_manual(name="Genotype",values = ann_colors$Genotype)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  labs(x="",y="Normalized Counts",
       title = paste("Gene with lowest adjusted p-value: ",g1,sep=""))+
  theme_bw()+
  theme(plot.title = element_text(family = "Arial", color="black",  size=20, hjust=0)) +
  theme(axis.title = element_text(family = "Arial", color="black", size=18))+
  theme(axis.text.x = element_text(family = "Arial",colour="black",size=14,angle = 32, hjust = 1))+
  theme(axis.text.y = element_text(family = "Arial",colour="black",size=14))+
  theme(legend.text = element_text(family = "Arial",colour="black",size=14))+
  theme(legend.title = element_blank())
#theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#    panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(paste("lowest_p.png",sep=""),path="./Figures/",
       width = 20, height = 15, units = "cm")
d<-plotCounts(ddsTxi, gene=ind2, intgroup="Genotype",
              main =g2,returnData = T)
p<-ggplot(d,aes(x=Genotype, y=count))
p+geom_boxplot(aes(fill=Genotype),alpha = 0.6)+
  geom_point(aes(fill=Genotype),colour="black",shape=21,
             position=position_jitter(w=0.05,h=0),
             na.rm=T,size=1.5)+
  guides(colour=guide_legend(ncol=2),fill=guide_legend(ncol=2))+
  scale_color_manual(name="Genotype",values = ann_colors$Genotype)+
  scale_fill_manual(name="Genotype",values = ann_colors$Genotype)+
  labs(x="",y="Normalized Counts",
       title = paste("Gene in with highest FC: ",g2))+
  theme_bw()+
  theme(plot.title = element_text(family = "Arial", color="black",  size=20, hjust=0)) +
  theme(axis.title = element_text(family = "Arial", color="black", size=18))+
  theme(axis.text.x = element_text(family = "Arial",colour="black",size=14,angle = 32, hjust = 1))+
  theme(axis.text.y = element_text(family = "Arial",colour="black",size=14))+
  theme(legend.text = element_text(family = "Arial",colour="black",size=14)) +
  theme(legend.title = element_blank())
#theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(paste("highest_FC.png",sep=""),path="./Figures",
       width = 20, height = 15, units = "cm")

#### Diagnostic Plots ####
# MA plot
png(filename = "./Figures/MAplot.png",width=15,height=10,units="cm",res=300)
plotMA(resWTvsKO, ylim=c(-4,4))
dev.off()

resultsNames(ddsTxi)
resLFC <- lfcShrink(ddsTxi, coef="Genotype_WT_vs_KO", type="apeglm")

png(filename = "./Figures/MAplot_LFC.png",width=15,height=10,units="cm",res=300)
plotMA(resLFC, ylim=c(-1,1))
dev.off()

#idx <- identify(resWTvsKO$baseMean, resWTvsKO$log2FoldChange)
#goi<-rownames(resWTvsKO)[idx]

mcols(resWTvsKO)$description
write.table(resWTvsKO[which(resWTvsKO$padj<0.05),],
            file="siggenes_WT_vs_KO.txt",col.names=T,row.names=F,quote=F,sep="\t")
write.table(resWTvsKO,file="allgenes_WT_vs_KO.txt",col.names=T,row.names=F,quote=F,sep="\t")

# resWTvsKO2<-resWTvsKO2[!is.na(resWTvsKO2$external_gene_name),]
# row.names(resWTvsKO2)<-resWTvsKO2$external_gene_name
# resultsReport <- HTMLReport(shortName ='WTvsKO_RNAseq_analysis_with_DESeq2',
#                          title ='RNA-seq analysis of differential expression using DESeq2',
#                          reportDirectory = "./reports")
# publish(resWTvsKO2,resultsReport, DataSet=ddsTxi, pvalueCutoff=0.01, lfc=0.5,
#         annotation.db=NULL, 
#         reportDir="./reports")
# finish(resultsReport)
# 

des2Report <- HTMLReport(shortName ='RNAseq_analysis_with_DESeq2',
                         title ='RNA-seq analysis of differential expression using DESeq2',
                         reportDirectory = "./reports")
publish(ddsTxi,resultsReport, n=1000, pvalueCutoff=0.05,
        annotation.db=NULL, factor = colData(ddsTxi)$Genotype,
        reportDir="./reports")
finish(des2Report)

#### Plot histogram of raw p-values in the results set
png(filename = "./Figures/hist_WTvsKO.png",width=15,height=10,units="cm",res=300)
hist(resWTvsKO$padj,xlab="Adjusted p-value",main="Adjusted p-value distribution")
dev.off()

cat("WT vs KO:\n",nrow(resWTvsKO[which(resWTvsKO$padj<0.05),]),
    nrow(resWTvsKO[which(resWTvsKO$padj<0.01),]),
    nrow(resWTvsKO[which(resWTvsKO$padj<0.002),]),"\n")

png(filename = "./Figures/hist_raw_WTvsKO.png",width=15,height=10,units="cm",res=300)
hist(resWTvsKO$pvalue,xlab="Raw p-value",main="Raw p-value distribution")
dev.off()

# with(resWTvsKO,plot(log2FoldChange,-log(padj),xlab="log2FoldChange",ylab="-log(p-value)",
#                     main="padj<0.002",pch=20,cex=0.5,ylim=c(0,30),col="#99999920"))
# with(resWTvsKO[which(resWTvsKO$padj<0.002),],points(log2FoldChange,-log(padj),
#                                                     col="#FF000040",pch=20,cex=0.5))
# 
# 
# resIDs<-sapply(strsplit(rownames(resWTvsKO),"[.]"),"[",1)
# genemap_idx2<-match(resIDs,genemap$ensembl_gene_id)
# results$external_gene_name<-genemap$external_gene_name[genemap_idx2]
# results$sig<-results$padj<0.05

#rownames(results)<-genemap$external_gene_name[genemap_idx]

results<-as.data.frame(resWTvsKO)
results<-results[!is.na(results$padj),]
results<-results[order(results$padj,results$log2FoldChange),]
sig<-results$padj<0.05
p <- ggplot(results, aes(log2FoldChange, -log10(padj))) +
  geom_point(ggplot2::aes(col = sig),alpha=0.5,show.legend = F) +
  scale_color_manual(values = c("black", "red")) +
  ggtitle("Volcano Plot") +
  geom_text_repel(data=results[1:10, ], aes(label=external_gene_name))+
  theme(plot.title = element_text(family = "Arial", 
                                  color="black",  size=20, hjust=0)) +
  theme(axis.title = element_text(family = "Arial", color="black", size=18))+
  theme(axis.text.x = element_text(family = "Arial",colour="black",size=14))+
  theme(axis.text.y = element_text(family = "Arial",colour="black",size=14))+
  theme(legend.text = element_text(family = "Arial",colour="black",size=14)) +
  theme(legend.title = element_blank())+
  theme(legend.text = element_text(family = "Arial"))+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))
p
ggsave("VolcanoPlot.png",path="./Figures",
       width = 20, height = 15, units = "cm")

results<-results[results$padj<0.01,]
results<-results[order(results$log2FoldChange,decreasing = T),]
goi<-row.names(results)[1:24]


stopifnot(all(goi %in% rownames(ddsTxi)))
tcounts <- t(assay(ddsTxi)[goi,]) %>%
  merge(colData(ddsTxi), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))
#tcounts %>% 
#  select(Row.names, Genotype, Batch, gene, expression) %>% 
#  head %>% 
#  knitr::kable()
tcountsIDs<-sapply(strsplit(tcounts$gene,"[.]"),"[",1)
genemap_idx3<-match(tcountsIDs,genemap$ensembl_gene_id)
tcounts$symbol<-genemap$external_gene_name[genemap_idx3]
tcounts$symbol[which(tcounts$symbol=="")]<-tcountsIDs[which(tcounts$symbol=="")]
p<- ggplot(tcounts, aes(Genotype, expression, fill=Genotype))
p+geom_boxplot() + 
  facet_wrap(~symbol, nrow=3) +
  geom_point(aes(fill=Genotype),colour="#666666",shape=21,
             position = position_dodge2(width=0.25,preserve = "total"))+
  #scale_y_log10()+
  scale_color_manual(name="Genotype",values = ann_colors$Genotype)+
  scale_fill_manual(name="Genotype",values = ann_colors$Genotype)+
  #guides(colour=guide_legend(ncol=2),fill=guide_legend(ncol=2))+
  xlab("")+
  labs(x="",
       y="Expression\n(Normalized Counts)", 
       fill="Genotype", 
       title="")+
  #theme_bw()+
  theme(axis.ticks.x = element_blank(),
        axis.ticks = element_line(colour="black"))+
  theme(plot.title = element_text(family = "Arial", color="black",  size=20, hjust=0)) +
  theme(axis.title = element_text(family = "Arial", color="black", size=18))+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(family = "Arial",colour="black",size=14))+
  theme(legend.text = element_text(family = "Arial",colour="black",size=14)) +
  theme(legend.title = element_blank())+
  theme(strip.text = element_text(family="Arial",color="black",size=14))+
  theme(legend.text = element_text(family = "Arial"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("TopGenes.png",path="./Figures",
       width = 20, height = 10, units = "cm")

write.csv(results,file="results_filtered.csv",row.names=T,quote=F)

##### Plots specific genes of interest #####
goi<-c("Ccne1","Dbf4","Mcm2","Mcm5","Cdc20","Cdk4","Cdc6","Cdk1","Mcm4",
       "Pkmyt1","Mcm6","Pcna","Rbl1","Ccne2","Mcm7","Mad2l1","Plk1","Chek1",
       "Ccnb2","Bub1b","Ccnb1","Mcm3","Cdc25c","Espl1")

goi_ensembl<-genemap$ensembl_gene_id[match(goi,genemap$external_gene_name)]

#rownames(ddsTxi)<-filteredIDs

tcounts<- t(assay(ddsTxi)[match(goi_ensembl,rownames(ddsTxi)),]/
              DESeq2::normalizationFactors(ddsTxi)[match(goi_ensembl,rownames(ddsTxi)),]) %>%
  merge(colData(ddsTxi), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi_ensembl)+1):ncol(.)) #11 genes

tcounts$symbol<-genemap$external_gene_name[match(tcounts$gene,genemap$ensembl_gene_id)]

p<- ggplot(tcounts, aes(Genotype, expression, fill=Genotype))
p+geom_boxplot() + 
  facet_wrap(~symbol, nrow=3) +
  geom_point(aes(fill=Genotype),colour="#666666",shape=21,
             position = position_dodge2(width=0.25,preserve = "total"))+
  #scale_y_log10()+
  scale_color_manual(name="Genotype",values = ann_colors$Genotype)+
  scale_fill_manual(name="Genotype",values = ann_colors$Genotype)+
  #guides(colour=guide_legend(ncol=2),fill=guide_legend(ncol=2))+
  xlab("")+
  labs(x="",
       y="Expression (Normalized Counts)", 
       fill="Condition", 
       title="")+
  #theme_bw()+
  theme(axis.ticks.x = element_blank(),
        axis.ticks = element_line(colour="black"))+
  theme(plot.title = element_text(family = "Arial", color="black",  size=20, hjust=0)) +
  theme(axis.title = element_text(family = "Arial", color="black", size=18))+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(family = "Arial",colour="black",size=14))+
  theme(legend.text = element_text(family = "Arial",colour="black",size=14)) +
  theme(legend.title = element_blank())+
  theme(strip.text = element_text(family="Arial",color="black",size=14))+
  theme(legend.text = element_text(family = "Arial"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("goi.png",path="./Figures",width = 20, height = 10, units = "cm")

#### Top 50 most different genes #####
select <- order(resWTvsKO$padj,decreasing = F)[1:50]

png(filename = "./Figures/QC/TopSig50Heatmap.png",width=15,height=20,units="cm",res=300)
pheatmap(assay(vsd.fast)[select,order(colData(vsd.fast)$Genotype)],
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "PiYG")))(100),
         cluster_rows=FALSE, show_rownames=T,
         cluster_cols=FALSE, annotation_col=df,
         annotation_colors = ann_colors[1])
dev.off()

g1<-resWTvsKO$external_gene_name[ind1]
ind2<-which.max(resWTvsKO$log2FoldChange)
g2<-resWTvsKO$external_gene_name[ind2]