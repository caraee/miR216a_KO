###### DESeq Analysis using tximport for scRNA-seq ######
#### Set up environment ####

library("BiocManager")
library("here")
library("biomaRt")
library("tximport")
library("DESeq2")
BiocManager::install("readr")

#### Load files ####
sampleTable<-data.frame(SampleID=c("SE38-13_S5", "SE38-16_S6","SE38-17_S7",
                                 "SE38-19_S8","SE38-20_S9","SE38-25_S10",
                                 "SE38-29_S11","SE38-30_S12"),
                        Genotype=c("WT","WT","KO","KO","KO","KO","WT","WT"))


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

#load("~/RProjects/ViaCyteRNASeq/txi.RData")


#Code used on cluster cluster
#library("tximport")
#library("DESeq2")
#library("BiocParallel")
#register(MulticoreParam(8))
#ddsTxiSC <- DESeq(ddsTxiSC,parallel=T,BPPARAM=MulticoreParam(8))
#save(txiSC,ddsTxiSC,file="txiSC.RData")

#possibly come back to this - maybe combine CellLine and CellStage?
sampleTable$Days<-gsub("S[4-7]D","",sampleTable$Stage)
sampleTable$Days<-gsub("HumanIslet","238",sampleTable$Days)
sampleTable$Days<-gsub("PEC01","238",sampleTable$Days)
sampleTable$Days<-as.numeric(sampleTable$Days)

sampleTable$CellStage<-gsub("D[0-9]+","",sampleTable$Stage)
sampleTable$CellStage<-gsub("HumanIslet","10",sampleTable$CellStage)

txiTC <- tximport(files[sampleTable$Condition!="HumanIslet"&sampleTable$Condition!="S"], type="salmon", tx2gene=tx2gene, ignoreTxVersion = T, #countsFromAbundance=c("lengthScaledTPM"),
                dropInfReps = T)

ddsTC <- DESeqDataSetFromTximport(txiTC,
                                   colData = sampleTable[sampleTable$Condition!="HumanIslet"&sampleTable$Condition!="S",],
                                   design = ~ CellLine+CellStage+Days+CellLine:Days:CellStage)


ddsTC<- DESeq(ddsTC)

library("scater")
browseVignettes("scater")
data("sc_example_counts")
data("sc_example_cell_info")
pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
rownames(pd) <- pd$Cell
example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
install.packages('Rtsne')
source('https://bioconductor.org/biocLite.R') 
biocLite('destiny')

pd<-new("AnnotatedDataFrame", data=sampleTable)
rownames(pd)<-pd$Sample
counts<-counts(ddsTxi)
counts<-counts[,row.names(pd)]
sceset<-newSCESet(countData=counts,phenoData = pd)
keep_feature <- rowSums(exprs(sceset) > 0) > 0 #30256 genes expressed in at least one cell
sceset <- sceset[keep_feature,] 
sceset <- calculateQCMetrics(sceset, feature_controls = 1:40)
scater_gui(sceset)

plotExpression(example_sceset, rownames(example_sceset)[7:12],
               x = "Mutation_Status", exprs_values = "counts", colour = "Cell_Cycle",
               show_median = TRUE, show_violin = FALSE,  xlab = "Mutation Status",
               log = TRUE)
plotFeatureData(example_sceset, aes(x = n_cells_exprs, y = pct_total_counts))
