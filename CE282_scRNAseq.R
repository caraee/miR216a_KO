##### Seurat #####
library("tidyverse")
library("here")
library("Seurat")
library("ggplot2")
library("cowplot")
library("RColorBrewer")

#### Prepare individual libraries #####
#### CE282_621665 #####
M621665 <- Read10X(data.dir = "~/RProjects/CE282_scRNAseq/data/CE282_621665/filtered_feature_bc_matrix")
M621665 <- CreateSeuratObject(counts = M621665, project = "CE282", min.cells = 3, min.features = 200)

M621665 #18445 features across 4954 samples
M621665[["Genotype"]]<-"KO"
M621665[["AnimalID"]]<-"M621665"
M621665[["percent.mt"]] <- PercentageFeatureSet(object = M621665, pattern = "^mt-")

png(filename = "./Figures/QC/M621665_pre_QC.png",width=15,height=15,units="cm",res=300)
  VlnPlot(M621665, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

summary(M621665$percent.mt)
summary(M621665$nFeature_RNA)
summary(M621665$nCount_RNA)

M621665 <- subset(x = M621665, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & 
                      percent.mt < 5)
M621665 #18445 features across 3122 samples - pretty good.

png(filename = "./Figures/QC/M621665_post_QC.png",width=15,height=15,units="cm",res=300)
  VlnPlot(M621665, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

save(M621665,file="./M621665.RData")
#rm(M621665)

#### CE282_621666 #####
M621666 <- Read10X(data.dir = "~/RProjects/CE282_scRNAseq/data/CE282_621666/filtered_feature_bc_matrix")
M621666 <- CreateSeuratObject(counts = M621666, project = "M621666", min.cells = 3, min.features = 200)
M621666 #18866 features across 11387 samples
M621666[["Genotype"]]<-"WT"
M621666[["AnimalID"]]<-"M621666"
M621666[["percent.mt"]] <- PercentageFeatureSet(object = M621666, pattern = "^mt-")

png(filename = "./Figures/QC/M621666_pre_QC.png",width=15,height=15,units="cm",res=300)
  VlnPlot(M621666, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

summary(M621666$percent.mt)
summary(M621666$nFeature_RNA)
summary(M621666$nCount_RNA)

M621666 <- subset(x = M621666, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & 
                    percent.mt < 5)
M621666 #18866 features across 5942 samples - excellent.

png(filename = "./Figures/QC/M621666_post_QC.png",width=15,height=15,units="cm",res=300)
  VlnPlot(M621666, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

save(M621666,file="./M621666.RData")

#### CE282_621669 #####
M621669 <- Read10X(data.dir = "~/RProjects/CE282_scRNAseq/data/CE282_621669/filtered_feature_bc_matrix")
M621669 <- CreateSeuratObject(counts = M621669, project = "M621669", min.cells = 3, min.features = 200)
M621669 #18740 features across 6908 samples
M621669[["Genotype"]]<-"KO"
M621669[["AnimalID"]]<-"M621669"
M621669[["percent.mt"]] <- PercentageFeatureSet(object = M621669, pattern = "^mt-")

png(filename = "./Figures/QC/M621669_pre_QC.png",width=15,height=15,units="cm",res=300)
  VlnPlot(M621669, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

summary(M621669$percent.mt)
summary(M621669$nFeature_RNA)
summary(M621669$nCount_RNA)

M621669 <- subset(x = M621669, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & 
                    percent.mt < 5)
M621669 #18740 features across 3668 samples.

png(filename = "./Figures/QC/M621669_post_QC.png",width=15,height=15,units="cm",res=300)
  VlnPlot(M621669, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

save(M621669,file="./M621669.RData")

#### CE282_621670 #####
M621670 <- Read10X(data.dir = "~/RProjects/CE282_scRNAseq/data/CE282_621670/filtered_feature_bc_matrix")
M621670 <- CreateSeuratObject(counts = M621670, project = "M621670", min.cells = 3, min.features = 200)
M621670 #18685 features across 5469 samples
M621670[["Genotype"]]<-"WT"
M621670[["AnimalID"]]<-"M621670"
M621670[["percent.mt"]] <- PercentageFeatureSet(object = M621670, pattern = "^mt-")

png(filename = "./Figures/QC/M621670_pre_QC.png",width=15,height=15,units="cm",res=300)
VlnPlot(M621670, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

summary(M621670$percent.mt)
summary(M621670$nFeature_RNA)
summary(M621670$nCount_RNA)

M621670 <- subset(x = M621670, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & 
                    percent.mt < 5)
M621670 #18685 features across 2302 samples

png(filename = "./Figures/QC/M621670_post_QC.png",width=15,height=15,units="cm",res=300)
VlnPlot(M621670, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

save(M621670,file="./M621670.RData")

##### Integrate datasets #####
load("./M621665.RData")
load("./M621666.RData")
load("./M621669.RData")
load("./M621670.RData")

options(future.globals.maxSize = 8000 * 1024^2)

islets.list <- merge(x=M621665,
                     y=c(M621666,M621669,M621670))
islets.list <- SplitObject(islets.list,split.by = "AnimalID")

for (i in 1:length(islets.list)) {
  islets.list[[i]] <- SCTransform(islets.list[[i]], verbose = FALSE)}

islets.features <- SelectIntegrationFeatures(object.list = islets.list, nfeatures = 3000)
islets.list <- PrepSCTIntegration(object.list = islets.list,
                                  anchor.features = islets.features,
                                  verbose = F)
islets.anchors <- FindIntegrationAnchors(object.list = islets.list,
                                         normalization.method = "SCT",
                                         anchor.features = islets.features,
                                         k.filter = 150,
                                         verbose = F)
islets.integrated <- IntegrateData(anchorset = islets.anchors, normalization.method = "SCT",
                                   verbose = F)

save(islets.integrated, file="islets.integrated.RData")

##### Full dataset QC #####
islets.integrated <- FindVariableFeatures(islets.integrated, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(islets.integrated), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(islets.integrated)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

islets.sce<-as.SingleCellExperiment(islets.integrated,assay="RNA")
gc()
library("scater")
islets.per.cell <- perCellQCMetrics(islets.sce, 
                                    subsets=list(Mito=grep("MT-", rownames(islets.sce))))
summary(islets.per.cell$sum)
summary(islets.per.cell$detected)
summary(islets.per.cell$subsets_Mito_percent)
colData(islets.sce) <- cbind(colData(islets.sce), islets.per.cell)

plotColData(islets.sce, x = "sum", y="detected", colour_by="Study") 
plotColData(islets.sce, x = "sum", y="detected", 
            other_fields="Study") + facet_wrap(~Study)


qc.stats <- quickPerCellQC(islets.per.cell, percent_subsets="subsets_Mito_percent")
colSums(as.matrix(qc.stats))
dim(filtered)

filtered <- islets.sce[,!qc.stats$discard] #19974 features across 15034 samples

png(filename = "./Figures/scater_highest_exprs.png",
    width=15,height=30,units="cm",res=300)
plotHighestExprs(islets.sce, exprs_values = "counts")
dev.off()

# keep_feature <- nexprs(islets.sce, byrow=TRUE) > 0
# islets.sce <- islets.sce[keep_feature,]
# dim(islets.sce)
# 
# example_sce <- logNormCounts(example_sce) # see below.
# vars <- getVarianceExplained(example_sce, 
#                              variables=c("tissue", "total mRNA mol", "sex", "age"))
# head(vars)
# plotExplanatoryVariables(vars)
# 
# example_sce <- logNormCounts(example_sce)
# assayNames(example_sce)
# 
# summary(librarySizeFactors(example_sce))
# 
# plotExpression(example_sce, rownames(example_sce)[1:6], x = "level1class")
# 
# plotExpression(example_sce, rownames(example_sce)[1:6],
#                x = rownames(example_sce)[10])
# 
# plotExpression(example_sce, rownames(example_sce)[1:6],
#                x = "level1class", colour_by="tissue")
# 
# plotExpression(example_sce, rownames(example_sce)[1:6])
# colData(example_sce) <- cbind(colData(example_sce), per.cell)

##### Dimensionality reduction #####
islets.integrated <- RunPCA(object = islets.integrated, verbose = FALSE)
islets.integrated <- RunUMAP(object = islets.integrated, dims = 1:30)
islets.integrated <- RunTSNE(object = islets.integrated, seed.use = 666,assay = "integrated",
                             check_duplicates = FALSE)

print(x = islets.integrated[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(islets.integrated, dims = c(1,4), reduction = "pca")
DimPlot(islets.integrated, reduction = "pca")
DimHeatmap(islets.integrated,  dims = 1:15, cells = 500, balanced = TRUE)

plot <- UMAPPlot(islets.integrated, group.by = c("AnimalID", "Genotype"))
png(filename = "./Figures/umap_mouse_genotype.png",
    width=35,height=17.5,units="cm",res=300)
plot & theme(legend.position = "top") & guides(color = guide_legend(nrow = 3, byrow = TRUE, 
                                                                     override.aes = list(size = 3)))
dev.off()

#plot most variable genes contributing to a PC for 500 most different cells - each column is a cell
png(filename = "./Figures/Diagnostics/PC_genes_human_islets.png",width=25,height=25,units="cm",res=300)
DimHeatmap(object = islets.integrated, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

png(filename = "./Figures/Diagnostics/PC_genes_islets_selectedPC.png",width=25,height=25,units="cm",res=300)
DimHeatmap(object = islets.integrated, dims = 1:6, cells = 500, balanced = TRUE)
dev.off()

save(islets.integrated,file="./islets_integrated.RData")

##### Cluster cells #####
DefaultAssay(islets.integrated) <- "integrated"
islets.integrated <- JackStraw(object = islets.integrated, num.replicate = 100)
islets.integrated <- ScoreJackStraw(object = islets.integrated, dims = 1:20)

png(filename = "./Figures/Diagnostics/jackstraw.png",width=25,height=15,units="cm",res=300)
JackStrawPlot(object = islets.integrated, dims = 1:20)
dev.off()

png(filename = "./Figures/Diagnostics/elbow.png",width=25,height=15,units="cm",res=300)
ElbowPlot(object = islets.integrated)
dev.off()

DefaultAssay(islets.integrated) <- "integrated"
islets.integrated <- FindNeighbors(islets.integrated, reduction = "pca", dims = 1:10, nn.eps = 0.5)
islets.integrated <- FindClusters(islets.integrated, resolution = 1.2, n.start = 23)
table(Idents(object = islets.integrated))

#save(islets.integrated,file="./islets_clustered.RData")
#load("./islets_clustered.RData")

cluster_cols2<-c("#343d00",
                 "#834efd",
                 "#d7ff44",
                 "#011b9e",
                 "#9acb00",
                 "#e953ff",
                 "#00df71",
                 "#ac00ba",
                 "#438500",
                 "#dc0086",
                 "#008e5a",
                 "#ff0072",
                 "#02edf7",
                 "#ff3401",
                 "#70a0ff",
                 "#ffdb44",
                 "#001756",
                 "#deffb0",
                 "#3b003a",
                 "#91ffe3",
                 "#790034",
                 "#a5dcff",
                 "#892f00",
                 "#a78cff",
                 "#827a00",
                 "#ff99ce",
                 "#002d00",
                 "#ff5f57",
                 "#029cb4",
                 "#ff9a62",
                 "#000f14",
                 "#ffecc1",
                 "#3a0019",
                 "#f9d2ff",
                 "#00536d")


png(filename = "./Figures/islets_clusters_umap.png",width=15,height=12,units="cm",res=300)
DimPlot(islets.integrated, reduction = "umap",cols = cluster_cols2,label = T,label.size = 5) + NoLegend()
dev.off()

plots <- UMAPPlot(islets.integrated, group.by = c("seurat_clusters", "Genotype"),
                  cols = cluster_cols2, label=T,label.size = 5,repel=T)

png(filename = "./Figures/umap_clusters_celltype.png",
    width=35,height=17.5,units="cm",res=300)
plots & theme(legend.position = "top") + NoLegend()
dev.off()

##### Normalize RNA for visualization/differential expression #####
DefaultAssay(islets.integrated) <- "RNA"
islets.integrated <- NormalizeData(islets.integrated)
islets.integrated <- FindVariableFeatures(islets.integrated, selection.method = "vst", nfeatures = 3000)
islets.integrated <- ScaleData(islets.integrated, 
                               vars.to.regress = "AnimalID",
                               features = rownames(islets.integrated))

sum(row.names(islets.integrated) %in% rownames(GetAssayData(islets.integrated, slot = 'scale.data')))

normalized_gene_expression<-islets.integrated@assays$RNA@data
colnames(normalized_gene_expression)<-paste(islets.integrated$AnimalID,colnames(normalized_gene_expression),sep="_")
write.csv(normalized_gene_expression,"scaled_gene_expression_matrix.csv",row.names = T)

###### Differential Expression #####
Idents(islets.integrated)<-"seurat_clusters"
# find markers for every cluster compared to all remaining cells, report only the positive ones
islets_markers <- FindAllMarkers(object = islets.integrated, only.pos = T, assay="RNA",
                                 slot = "data",
                                 min.pct = 0.05, logfc.threshold = 0.15,
                                 return.thresh=0.05)

islets_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10_islets
islets_markers %>% group_by(cluster) -> islets_cluster_identifiers

write.csv(islets_cluster_identifiers,"islets_cluster_identifiers.csv",row.names = F)

png(filename = "./Figures/islets_heatmap_top10_clusters.png",width=25,height=40,units="cm",res=300)
DoHeatmap(object = islets.integrated, group.colors = cluster_cols2, 
          features = top10_islets$gene, label=F,
          assay="RNA") #+NoLegend()
dev.off()

# are some clusters likely doublets?
png(filename = "./Figures/clusters_QC.png",width=35,height=20,units="cm",res=300)
VlnPlot(islets.integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

#### Rename clusters based on identifying genes ####
Idents(islets.integrated) <- "seurat_clusters"
#Significantly upregulated genes
#α - Gcg
#β - Ins1, Ins2, and Iapp
#δ - Sst
#ε - Ghrl
#PPY - Ppy
#Ductal - Krt19
#Acinar - Cpa1
#Endothelial - Pecam1
#Macrophage - Cd68
#Stellate - Rgs5 (quiescent), Pdgfrb (activated)

new.cluster.ids <- c("β.2","β.1","β.2","β.2","4","α/PPY","β.1", "β.1","8",
                     "endothelial","β/acinar","δ/PPY","acinar","δ/PPY",
                     "β/acinar","macrophage","acinar","macrophage","18","PPY",
                     "acinar","acinar","α/PPY","stellate","24","25",
                     "endothelial", "α/δ/PPY","stellate","ductal")
names(new.cluster.ids) <- levels(islets.integrated$seurat_clusters)
islets.integrated <- RenameIdents(islets.integrated, new.cluster.ids)
islets.integrated$celltype <- Idents(islets.integrated)

islets.integrated$celltype.genotype <- paste(Idents(islets.integrated), islets.integrated$Genotype, sep = "_")
Idents(islets.integrated) <- "celltype"

#levels(islets.integrated)<-c("0","2","3","8","Ductal-like","β-like","4","α-like","γ/ε-like")

png(filename = "./Figures/celltype_umap_labelled.png",width=15,height=12,units="cm",res=300)
DimPlot(islets.integrated, cols=cluster_cols2, 
        reduction = "umap",
        label.size = 5,
        label=T,repel = T)+NoLegend()
dev.off()

##### Propotions of cells in each celltype ####
table(islets.integrated$Genotype)
table(islets.integrated$celltype.genotype)

c<-which(islets.integrated$seurat_clusters%in%"7")
summary(islets.integrated@assays$RNA@data["Ins1",c]) #Ins1>6.986
summary(islets.integrated@assays$RNA@data["Ins2",c]) #Ins2>8.223
summary(islets.integrated@assays$RNA@data["Iapp",c]) #Iapp>4.512

DefaultAssay(islets.integrated) <- "RNA"
Idents(islets.integrated) <- "seurat_clusters"
beta_thresh<-WhichCells(object=islets.integrated, expression=Ins1>6.986&Ins2>8.223&Iapp>4.512)
beta_clust<-WhichCells(object=islets.integrated, idents = c("0","1","2","3","4","6",
                                                            "7","8","10","14","18","25"))
betas<-intersect(beta_thresh,beta_clust)

islets.integrated@assays$RNA@counts["GRCh38.98---LEP",lep]

Idents(islets.integrated)<-"celltype"
png(filename = "./Figures/human_islets/islets_umap_lep.png",width=15,height=10,units="cm",res=300)
DimPlot(islets.integrated, cells.highlight = lep, 
        reduction = "umap",label=T,repel=T)+NoLegend()
dev.off()

png(filename = "./Figures/human_islets/islets_feature_umap_lep_lepr.png",width=30,height=15,units="cm",res=300)
FeaturePlot(object = islets.integrated, reduction = "umap",
            features = c("GRCh38.98---LEP","GRCh38.98---LEPR"), pt.size = 0.2)
dev.off()

png(filename = "./Figures/human_islets/islets_feature_umap_lep.png",width=15,height=10,units="cm",res=300)
FeaturePlot(object = islets.integrated, reduction = "umap",
            features = c("GRCh38.98---PRLR"), pt.size = 0.2)
dev.off()


##### Ridge Plots for Genes of Interest #####
Idents(islets.integrated)<-"seurat_clusters"
other<-c("Krt19","Cpa1","Pecam1","Krt19","Ptf1a","Sparc","Pecam1","Rgs5","Pdgfrb","Sox10","Ptprc",
         "Lyz2","Flt3","Trac","Iglc3")

DefaultAssay(islets.integrated) <- "RNA"
png(filename = "./Figures/goi_UMAP_other1.png",width=25,height=21,units="cm",res=300)
FeaturePlot(object = islets.integrated, 
            order = T,
            features = other[1:4], pt.size = 0.2) &
  theme( plot.title = element_text( face = "italic") )
dev.off()


png(filename = "./Figures/other_ridge1.png",width=30,height=15,units="cm",res=300)
RidgePlot(object = islets.integrated, cols=cluster_cols2, #group.by ="Genotype",
          features = other[1:3], log = TRUE,sort = "increasing")
dev.off()

png(filename = "./Figures/other_ridge2.png",width=30,height=15,units="cm",res=300)
RidgePlot(object = islets.integrated, cols=cluster_cols2,#group.by ="Genotype", 
          features = other[4:6], log = TRUE,sort = "increasing")
dev.off()

png(filename = "./Figures/other_ridge3.png",width=30,height=15,units="cm",res=300)
RidgePlot(object = islets.integrated, cols=cluster_cols2,#group.by ="Genotype", 
          features = other[7:9], log = TRUE,sort = "increasing")
dev.off()

png(filename = "./Figures/other_ridge4.png",width=30,height=15,units="cm",res=300)
RidgePlot(object = islets.integrated, cols=cluster_cols2,#group.by ="Genotype", 
          features = other[10:12], log = TRUE,sort = "increasing")
dev.off()

#### Enterochromaffin #####
ECL<-c("Ddc","Tph1","Lmx1a","Slc18a1","Tac1","Adra2a","Fev","Cxcl14",
       "Nkx6-1","Ins1","Ins2","G6pc2","Nptx2","Isl1","Pdx1")

png(filename = "./Figures/goi_UMAP_ECL1.png",width=20,height=28,units="cm",res=300)
FeaturePlot(object = islets.integrated, features = ECL[1:5], pt.size = 0.2, order = T)
dev.off()

png(filename = "./Figures/goi_UMAP_ECL2.png",width=20,height=28,units="cm",res=300)
FeaturePlot(object = islets.integrated, features = ECL[6:11], pt.size = 0.2,order=T)
dev.off()

png(filename = "./Figures/goi_UMAP_ECL3.png",width=20,height=21,units="cm",res=300)
FeaturePlot(object = islets.integrated, features = ECL[12:15], pt.size = 0.2,order=T)
dev.off()

##### Endocrine markers #####
endocrine<-c("Ins1","Ins2","Iapp","Gcg","Hhex","Sst","Ppy","Ghrl")

sst<-c("SST","LAMC2", "MAP3K15","SLC7A7", "TBX2B")
sst<-paste("",sst,sep="")

DefaultAssay(islets.integrated) <- "RNA"
png(filename = "./Figures/goi_UMAP_endocrine1.png",width=25,height=21,units="cm",res=300)
FeaturePlot(object = islets.integrated, 
            order = T,
            features = endocrine[1:4], pt.size = 0.2) &
  theme( plot.title = element_text( face = "italic") )
dev.off()

png(filename = "./Figures/goi_UMAP_endocrine2.png",width=25,height=21,units="cm",res=300)
FeaturePlot(object = islets.integrated, 
            order = T,
            features = endocrine[5:8], pt.size = 0.2) &
  theme( plot.title = element_text( face = "italic") )
dev.off()

png(filename = "./Figures/goi_UMAP_beta_acinar.png",width=25,height=21,units="cm",res=300)
FeaturePlot(object = islets.integrated, 
            order = T,
            features = c("Ins1","Ins2","Iapp","Cpa1"), pt.size = 0.2) &
  theme( plot.title = element_text( face = "italic") )
dev.off()

#contamination from Ins1, Ins2 and Iapp is highly prevalent

png(filename = "./Figures/goi_UMAP_endocrine2.png",width=20,height=14,units="cm",res=300)
FeaturePlot(object = islets.integrated, split.by="Genotype",features = endocrine[c(4:6)], pt.size = 0.2)
dev.off()

png(filename = "./Figures/endocrine_ridge1.png",width=30,height=15,units="cm",res=300)
RidgePlot(object = islets.integrated, cols=cluster_cols2, #group.by ="Genotype",
          features = endocrine[1:3], log = TRUE,sort = "increasing")
dev.off()

png(filename = "./Figures/endocrine_ridge2.png",width=20,height=15,units="cm",res=300)
RidgePlot(object = islets.integrated, cols=cluster_cols2,#group.by ="Genotype", 
          features = endocrine[4:6], log = TRUE,sort = "increasing")
dev.off()

png(filename = "./Figures/endocrine_ridge3.png",width=20,height=15,units="cm",res=300)
RidgePlot(object = islets.integrated, cols=cluster_cols2,#group.by ="Genotype", 
          features = endocrine[6:7], log = TRUE,sort = "increasing")
dev.off()

png(filename = "./Figures/endocrine_splitdot1.png",width=20,height=10,units="cm",res=300)
DotPlot(object = islets.integrated, group.by="Genotype",features = endocrine[1:3]) + RotatedAxis()
dev.off()

png(filename = "./Figures/endocrine_splitdot2.png",width=20,height=10,units="cm",res=300)
DotPlot(object = islets.integrated, group.by="Genotype",features = endocrine[c(4,5)]) + RotatedAxis()
dev.off()


##### Stromal cells #####
stromal<-c("Hic1","Pdgfra","Tgfb1","Col1a1","Vim","Ccn2","Fn1")

png(filename = "./Figures/goi_UMAP_stromal1.png",width=20,height=28,units="cm",res=300)
FeaturePlot(object = islets.integrated, split.by="Genotype",features = stromal[1:4], pt.size = 0.2)
dev.off()

png(filename = "./Figures/goi_UMAP_stromal2.png",width=20,height=21,units="cm",res=300)
FeaturePlot(object = islets.integrated, split.by="Genotype",features = stromal[5:7], pt.size = 0.2)
dev.off()

png(filename = "./Figures/stromal_ridge.png",width=30,height=15,units="cm",res=300)
RidgePlot(object = islets.integrated, cols=cluster_cols, idents = levels(islets.integrated)[1:10],#group.by ="Genotype", 
          features = stromal[1:6], slot = 'counts', log = TRUE,
          pt.size = 0)
dev.off()

##### Immune cells #####
immune<-c("CD4","CD8A","PTPRC","CD68","CD163A","CD14","ITGAM","ITGAX")
immune<-paste("",immune,sep="")

png(filename = "./Figures/immune_ridge.png",width=30,height=15,units="cm",res=300)
RidgePlot(object = islets.integrated, cols=cluster_cols, idents = levels(islets.integrated)[1:10],#group.by ="Genotype", 
          features = immune[1:8], slot = 'counts', log = TRUE,
          pt.size = 0)
dev.off()

##### Insulin biosynthesis #####
biosynthesis<-c("Pcsk1","Pcsk2",
                "Abcc8","Kcnj11","Pclo",
                "Rims2","Sytl4","Vamp2",
                "Gjd2","Slc30a8")

png(filename = "./Figures/goi_UMAP_GLP1.png",width=20,height=21,units="cm",res=300)
FeaturePlot(object = islets.integrated, split.by="Genotype",features = c(biosynthesis[1:2],endocrine[3]), pt.size = 0.2)
dev.off()

png(filename = "./Figures/goi_UMAP_biosynthesis1.png",width=20,height=21,units="cm",res=300)
FeaturePlot(object = islets.integrated, split.by="Genotype",features = biosynthesis[1:2], pt.size = 0.2)
dev.off()

png(filename = "./Figures/biosynthesis_ridge1.png",width=30,height=15,units="cm",res=300)
RidgePlot(object = islets.integrated, cols=cluster_cols2,#group.by ="Genotype", 
          features = biosynthesis[1:3], log = TRUE,sort="increasing")
dev.off()

png(filename = "./Figures/goi_UMAP_biosynthesis2.png",width=20,height=28,units="cm",res=300)
FeaturePlot(object = islets.integrated, split.by="Genotype",features = biosynthesis[4:6], pt.size = 0.2)
dev.off()

png(filename = "./Figures/biosynthesis_ridge2.png",width=30,height=15,units="cm",res=300)
RidgePlot(object = islets.integrated, cols=cluster_cols,#group.by ="Genotype", 
          features = biosynthesis[4:6], slot = 'counts', log = TRUE)
dev.off()

png(filename = "./Figures/goi_UMAP_biosynthesis3.png",width=20,height=21,units="cm",res=300)
FeaturePlot(object = islets.integrated, split.by="Genotype",features = biosynthesis[7:9], pt.size = 0.2)
dev.off()

png(filename = "./Figures/biosynthesis_ridge3.png",width=30,height=15,units="cm",res=300)
RidgePlot(object = islets.integrated, cols=cluster_cols,#group.by ="Genotype", 
          features = biosynthesis[7:9], slot = 'counts', log = TRUE)
dev.off()

png(filename = "./Figures/goi_UMAP_biosynthesis4.png",width=20,height=7,units="cm",res=300)
FeaturePlot(object = islets.integrated, split.by="Genotype",features = biosynthesis[10], pt.size = 0.2)
dev.off()

png(filename = "./Figures/biosynthesis_ridge4.png",width=15,height=15,units="cm",res=300)
RidgePlot(object = islets.integrated, cols=cluster_cols,#group.by ="Genotype", 
          features = biosynthesis[10], slot = 'counts', log = TRUE)+NoLegend()
dev.off()


##### Endothelial cells #####
#from https://www.physiology.org/doi/pdf/10.1152/physiolgenomics.00186.2002
endothelial<-c("Pecam1","Cdh5","C1QBP","MMRN1","EFEMP1",
               "VWF","ESM1","EDG1","NTN4") #no VWF, EDG1
endothelial<-paste("",endothelial,sep="")

png(filename = "./Figures/goi_UMAP_endothelial1.png",width=20,height=21,units="cm",res=300)
FeaturePlot(object = islets.integrated, split.by="Genotype",features = endothelial[1:3], pt.size = 0.2)
dev.off()

png(filename = "./Figures/goi_UMAP_endothelial2.png",width=20,height=21,units="cm",res=300)
FeaturePlot(object = islets.integrated, split.by="Genotype",features = endothelial[4:6], pt.size = 0.2)
dev.off()

png(filename = "./Figures/goi_UMAP_endothelial3.png",width=20,height=21,units="cm",res=300)
FeaturePlot(object = islets.integrated, split.by="Genotype",features = endothelial[7:9], pt.size = 0.2)
dev.off()


##### Metabolism #####
metab<-c("SLC2A1","GCK","HK1","ALDOB","GAPDH","ENO1",
         "LDHA","PC","COX6A1") #no ATP5A1
metab<-paste("",metab,sep="")

png(filename = "./Figures/goi_UMAP_metab1.png",width=20,height=21,units="cm",res=300)
FeaturePlot(object = islets.integrated, split.by="Genotype",features = metab[1:3], pt.size = 0.2)
dev.off()

png(filename = "./Figures/metab_ridge1.png",width=30,height=15,units="cm",res=300)
RidgePlot(object = islets.integrated, cols=cluster_cols,#group.by ="Genotype", 
          features = metab[1:3], slot = 'counts', log = TRUE,
          pt.size = 0)
dev.off()

png(filename = "./Figures/goi_UMAP_metab2.png",width=20,height=21,units="cm",res=300)
FeaturePlot(object = islets.integrated, split.by="Genotype",features = metab[4:6], pt.size = 0.2)
dev.off()

png(filename = "./Figures/metab_ridge2.png",width=30,height=15,units="cm",res=300)
RidgePlot(object = islets.integrated, cols=cluster_cols,#group.by ="Genotype", 
          features = metab[4:6], slot = 'counts', log = TRUE,
          pt.size = 0)
dev.off()

png(filename = "./Figures/goi_UMAP_metab3.png",width=20,height=21,units="cm",res=300)
FeaturePlot(object = islets.integrated, split.by="Genotype",features = metab[7:9], pt.size = 0.2)
dev.off()

png(filename = "./Figures/metab_ridge3.png",width=30,height=15,units="cm",res=300)
RidgePlot(object = islets.integrated, cols=cluster_cols,#group.by ="Genotype", 
          features = metab[7:9], slot = 'counts', log = TRUE,
          pt.size = 0)
dev.off()

##### Beta Cell Maturation Markers #####
mature<-c("Neurod1","Neurog3","Mafa","Mafb","Slc16a1","Ucn3",
          "Dlk1","Npy","Oat","Pdgfra","Igfbp4","Mylk","Cfap12")

png(filename = "./Figures/goi_UMAP_mature1.png",width=20,height=21,units="cm",res=300)
FeaturePlot(object = islets.integrated, split.by="Genotype",features = mature[1:3], pt.size = 0.2)
dev.off()

png(filename = "./Figures/mature_ridge1.png",width=30,height=15,units="cm",res=300)
RidgePlot(object = islets.integrated, cols=cluster_cols,#group.by ="Genotype", 
          features = mature[1:3], slot = 'counts', log = TRUE,
          pt.size = 0)
dev.off()

png(filename = "./Figures/goi_UMAP_mature2.png",width=20,height=21,units="cm",res=300)
FeaturePlot(object = islets.integrated, split.by="Genotype",features = mature[4:6], pt.size = 0.2)
dev.off()

png(filename = "./Figures/goi_UMAP_mature3.png",width=20,height=21,units="cm",res=300)
FeaturePlot(object = islets.integrated, split.by="Genotype",features = mature[7:9], pt.size = 0.2)
dev.off()

png(filename = "./Figures/goi_UMAP_mature4.png",width=20,height=21,units="cm",res=300)
FeaturePlot(object = islets.integrated, split.by="Genotype",features = mature[10:12], pt.size = 0.2)
dev.off()

##### Fat #####
panfat<-c("FABP4", "PPARA", "GPDH", "SCD1")
panfat<-paste("",panfat,sep="")
adipose_expansion<-c("MEST", "SFRP5", "CAV1", "BMP3")
adipose_expansion<-paste("",adipose_expansion,sep="")
brown<-c("UCP1", "PGC1A", "PPARA") 
brown<-paste("",brown,sep="")
lipid_homeostasis<-c("ACLY", "ACACA", "SCL25A", "ELOVL6") 
lipid_homeostasis<-paste("",lipid_homeostasis,sep="")

FeaturePlot(object = islets.integrated, split.by="Genotype",
            slot="data",
            features = "Lep", pt.size = 0.2)

png(filename = "./Figures/goi_UMAP_panfat.png",width=20,height=21,units="cm",res=300)
FeaturePlot(object = islets.integrated, split.by="Genotype",features = panfat, pt.size = 0.2)
dev.off()

png(filename = "./Figures/ridgeplot_panfat.png",width=30,height=15,units="cm",res=300)
RidgePlot(object = islets.integrated, cols=cluster_cols,group.by ="celltype", 
          features = panfat, slot = 'counts', log = TRUE,
          pt.size = 0)
dev.off()

png(filename = "./Figures/splitdot_panfat.png",width=20,height=10,units="cm",res=300)
DotPlot(object = islets.integrated, group.by="celltype",features = panfat) + RotatedAxis()
dev.off()


png(filename = "./Figures/goi_UMAP_adipose_expansion.png",width=20,height=21,units="cm",res=300)
FeaturePlot(object = islets.integrated, split.by="Genotype",features = adipose_expansion, pt.size = 0.2)
dev.off()

png(filename = "./Figures/ridgeplot_adipose_expansion.png",width=30,height=15,units="cm",res=300)
RidgePlot(object = islets.integrated, cols=cluster_cols,group.by ="celltype", 
          features = adipose_expansion, slot = 'counts', log = TRUE,
          pt.size = 0)
dev.off()

png(filename = "./Figures/splitdot_adipose_expansion.png",width=20,height=10,units="cm",res=300)
DotPlot(object = islets.integrated, group.by="celltype",features = adipose_expansion) + RotatedAxis()
dev.off()

#Not useful - only PPARA detected
png(filename = "./Figures/goi_UMAP_brown.png",width=20,height=21,units="cm",res=300)
FeaturePlot(object = islets.integrated, split.by="Genotype",features = brown, pt.size = 0.2)
dev.off()

png(filename = "./Figures/ridgeplot_brown.png",width=30,height=15,units="cm",res=300)
RidgePlot(object = islets.integrated, cols=cluster_cols,group.by ="celltype", 
          features = brown, slot = 'counts', log = TRUE,
          pt.size = 0)
dev.off()

png(filename = "./Figures/splitdot_brown.png",width=20,height=10,units="cm",res=300)
DotPlot(object = islets.integrated, group.by="celltype",features = brown) + RotatedAxis()
dev.off()


png(filename = "./Figures/goi_UMAP_lipid_homeostasis.png",width=20,height=21,units="cm",res=300)
FeaturePlot(object = islets.integrated, split.by="Genotype",features = lipid_homeostasis, pt.size = 0.2)
dev.off()

png(filename = "./Figures/ridgeplot_lipid_homeostasis.png",width=30,height=15,units="cm",res=300)
RidgePlot(object = islets.integrated, cols=cluster_cols,group.by ="celltype", 
          features = lipid_homeostasis, slot = 'counts', log = TRUE,
          pt.size = 0)
dev.off()

png(filename = "./Figures/splitdot_lipid_homeostasis.png",width=20,height=10,units="cm",res=300)
DotPlot(object = islets.integrated, group.by="celltype",features = lipid_homeostasis) + RotatedAxis()
dev.off()


#### Differential expression by genotype and seurat cluster #####
Idents(islets.integrated) <- "celltype"
png(filename = "./Figures/celltype_Genotype_umap_labelled.png",width=30,height=12,units="cm",res=300)
DimPlot(islets.integrated, cols=cluster_cols2, split.by="Genotype", 
        reduction = "umap",
        label.size = 5,
        label=T,repel=T)+NoLegend()
dev.off()

#theme_set(theme_cowplot())

Idents(islets.integrated) <- "seurat_clusters"
png(filename = "./Figures/clusters_Genotype_umap_labelled.png",width=30,height=12,units="cm",res=300)
DimPlot(islets.integrated, cols=cluster_cols2, split.by="Genotype", 
        reduction = "umap",
        label.size = 5,
        label=T,repel=T)+NoLegend()
dev.off()


islets.integrated$cluster.genotype <- paste(Idents(islets.integrated), islets.integrated$Genotype, sep = "_")
Idents(islets.integrated) <- "cluster.genotype"
islets.integrated$celltype.donor <- paste(islets.integrated$celltype, islets.integrated$AnimalID, sep = "_")

levels(islets.integrated)<-levels(islets.integrated)[order(levels(islets.integrated))]
table(Idents(islets.integrated))

Idents(islets.integrated) <- "cluster.genotype"
for (i in levels(islets.integrated$seurat_clusters)) {
  ident.1 = paste(i,"KO",sep="_")
  ident.2 = paste(i,"WT",sep="_")
  cell.diffs <- FindMarkers(islets.integrated,
                              assay="RNA",
                              slot = "data",
                              min.pct = 0.1, logfc.threshold = 0.25,
                              return.thresh=0.1,
                              ident.1 = ident.1, 
                              ident.2 = ident.2, 
                              verbose = FALSE)
  cell.diffs<-cell.diffs[cell.diffs$p_val_adj<0.05,]
  cell.diffs<-cell.diffs[order(cell.diffs$avg_logFC,cell.diffs$p_val_adj),]
  # head(cell.diffs)
  # tail(cell.diffs)
  write.csv(cell.diffs,paste("cluster",i,"identifiers.csv",sep="_"),row.names = T)
  
  c<-which(Idents(islets.integrated)%in%c(ident.1,ident.2))
  if (length(cell.diffs$p_val_adj)>50) {feat<-c(head(row.names(cell.diffs),25),tail(row.names(cell.diffs),25))} else {
    feat<-row.names(cell.diffs)
  }
  
  if (length(cell.diffs$p_val_adj)>0) {
  sub<-subset(x = islets.integrated, cells = c, features=feat, idents=c(ident.1,ident.2))
  
  png(filename = paste("./Figures/heatmap_cluster",i,"cells.png",sep="_"),width=15,height=15,units="cm",res=300)
  DoHeatmap(object = sub, group.colors = cluster_cols2, features = feat, 
            slot="data",
            assay = "RNA",
            angle=0,label = F)
  dev.off()
  }
}

#why isn't the heatmap saving?

#### Differential expression by celltype and cluster #####
islets.integrated$celltype.genotype <- paste(Idents(islets.integrated), islets.integrated$Genotype, sep = "_")
Idents(islets.integrated) <- "celltype.genotype"

ident.1 = "β.1_KO"
ident.2 = "β.1_WT"
b1.cell.diffs <- FindMarkers(islets.integrated,
                          assay="RNA",
                          slot = "data",
                          min.pct = 0.1, logfc.threshold = 0.25,
                          return.thresh=0.1,
                          ident.1 = ident.1, 
                          ident.2 = ident.2, 
                          verbose = FALSE)
b1.cell.diffs<-b1.cell.diffs[b1.cell.diffs$p_val_adj<0.05,]
b1.cell.diffs<-b1.cell.diffs[order(b1.cell.diffs$avg_logFC,b1.cell.diffs$p_val_adj),]
# head(cell.diffs)
# tail(cell.diffs)
write.csv(b1.cell.diffs,"b1_identifiers.csv",row.names = T)

c<-which(Idents(islets.integrated)%in%c(ident.1,ident.2))
if (length(b1.cell.diffs$p_val_adj)>50) {feat<-c(head(row.names(b1.cell.diffs),25),tail(row.names(b1.cell.diffs),25))} else {
  feat<-row.names(b1.cell.diffs)
}
sub<-subset(x = islets.integrated, cells = c, features=feat, idents=c(ident.1,ident.2))

png(filename ="./Figures/heatmap_cluster_b1_cells.png",width=15,height=10,units="cm",res=300)
DoHeatmap(object = sub, group.colors = c("#ff8000",
                                         "#000080"), features = feat, 
            slot="data",
            assay = "RNA",
            angle=0,label = F)
dev.off()

ident.1 = "β.2_KO"
ident.2 = "β.2_WT"
b2.cell.diffs <- FindMarkers(islets.integrated,
                             assay="RNA",
                             slot = "data",
                             min.pct = 0.1, logfc.threshold = 0.25,
                             return.thresh=0.1,
                             ident.1 = ident.1, 
                             ident.2 = ident.2, 
                             verbose = FALSE)
b2.cell.diffs<-b2.cell.diffs[b2.cell.diffs$p_val_adj<0.05,]
b2.cell.diffs<-b2.cell.diffs[order(b2.cell.diffs$avg_logFC,b2.cell.diffs$p_val_adj),]
# head(cell.diffs)
# tail(cell.diffs)
write.csv(b2.cell.diffs,"b2_identifiers.csv",row.names = T)

c<-which(Idents(islets.integrated)%in%c(ident.1,ident.2))
if (length(b2.cell.diffs$p_val_adj)>50) {feat<-c(head(row.names(b2.cell.diffs),25),tail(row.names(b2.cell.diffs),25))} else {
  feat<-row.names(b2.cell.diffs)
}
sub<-subset(x = islets.integrated, cells = c, features=feat, idents=c(ident.1,ident.2))

png(filename ="./Figures/heatmap_cluster_b2_cells.png",width=15,height=15,units="cm",res=300)
DoHeatmap(object = sub, group.colors = c("#ff8000",
                                         "#000080"), features = feat, 
          slot="data",
          assay = "RNA",
          angle=0,label = F)
dev.off()


Idents(islets.integrated) <- "celltype.donor"
c<-which(Idents(islets.integrated)%in%c("β.2_M621665","β.2_M621666",
                                        "β.2_M621669","β.2_M621670"))
feat<-c("Ins1","Ins2","Iapp",row.names(b2.cell.diffs))
sub<-subset(x = islets.integrated, cells = c, features=feat, idents=c("β.2_M621665","β.2_M621666",
                                                                      "β.2_M621669","β.2_M621670"))
levels(sub)<-levels(sub)[c(1,3,2,4)]

png(filename = "./Figures/heatmap_betacells_donor.png",width=15,height=15,units="cm",res=300)
DoHeatmap(object = sub, features=feat, 
          slot="data",
          label=F, assay="RNA",angle=0)
dev.off()

Idents(islets.integrated) <- "celltype"
ident.1 = "acinar_KO"
ident.2 = "acinar_WT"
acinar.cell.diffs <- FindMarkers(islets.integrated,
                             assay="RNA",
                             slot = "data",
                             min.pct = 0.1, logfc.threshold = 0.25,
                             return.thresh=0.1,
                             ident.1 = ident.1, 
                             ident.2 = ident.2, 
                             verbose = FALSE)
acinar.cell.diffs<-acinar.cell.diffs[acinar.cell.diffs$p_val_adj<0.05,]
acinar.cell.diffs<-acinar.cell.diffs[order(acinar.cell.diffs$avg_logFC,acinar.cell.diffs$p_val_adj),]
# head(cell.diffs)
# tail(cell.diffs)
write.csv(acinar.cell.diffs,"acinar_identifiers.csv",row.names = T)

c<-which(Idents(islets.integrated)%in%c(ident.1,ident.2))
if (length(acinar.cell.diffs$p_val_adj)>50) {feat<-c(head(row.names(acinar.cell.diffs),25),tail(row.names(acinar.cell.diffs),25))} else {
  feat<-row.names(acinar.cell.diffs)
}
sub<-subset(x = islets.integrated, cells = c, features=feat, idents=c(ident.1,ident.2))

png(filename ="./Figures/heatmap_cluster_acinar_cells.png",width=15,height=15,units="cm",res=300)
DoHeatmap(object = sub, group.colors = c("#ff8000",
                                         "#000080"), features = feat, 
          slot="data",
          assay = "RNA",
          angle=0,label = F)
dev.off()


Idents(islets.integrated) <- "celltype.donor"
c<-which(Idents(islets.integrated)%in%c("acinar_M621665","acinar_M621666",
                                        "acinar_M621669","acinar_M621670"))
sub<-subset(x = islets.integrated, cells = c, features=feat, idents=c("acinar_M621665","acinar_M621666",
                                                                      "acinar_M621669","acinar_M621670"))
levels(sub)<-levels(sub)[c(1,3,2,4)]

png(filename = "./Figures/heatmap_acinarcells_donor.png",width=15,height=15,units="cm",res=300)
DoHeatmap(object = sub, features=feat, 
          slot="data",
          label=F, assay="RNA",angle=0)
dev.off()


#### Look at genes of interest by genotype and cluster #####
png(filename = "./Figures/rat_celltype_Genotype_clusters_labelled.png",width=25,height=15,units="cm",res=300)
DimPlot(islets.integrated, reduction = "umap", label = T,
        repel=T,cols = cluster_cols, pt.size = 0.5)#+NoLegend()
dev.off()

names(cluster_cols)<-levels(islets.integrated)

png(filename = "./Figures/endocrine_celltype_Genotype_ridge1.png",width=30,height=10,units="cm",res=300)
RidgePlot(object = islets.integrated, cols=cluster_cols, idents=levels(islets.integrated)[13:17],#group.by ="Genotype",
          features = endocrine[1:3], slot = 'counts', log = TRUE,
          pt.size = 0)
dev.off()

png(filename = "./Figures/endocrine_celltype_Genotype_ridge2.png",width=20,height=10,units="cm",res=300)
RidgePlot(object = islets.integrated, cols=cluster_cols,idents=levels(islets.integrated)[13:17],#group.by ="Genotype", 
          features = endocrine[4:5], slot = 'counts', log = TRUE,
          pt.size = 0)
dev.off()


png(filename = "./Figures/stromal_celltype_Genotype_ridge.png",width=20,height=5,units="cm",res=300)
RidgePlot(object = islets.integrated, idents=levels(islets.integrated)[1:6],cols=cluster_cols,#group.by ="Genotype", 
          features = stromal, slot = 'counts', log = F,
          pt.size = 0)
dev.off()

png(filename = "./Figures/mature_celltype_Genotype_ridge1.png",width=30,height=15,units="cm",res=300)
RidgePlot(object = islets.integrated, cols=cluster_cols,#group.by ="Genotype", 
          features = mature[1:3], slot = 'counts', log = TRUE,
          pt.size = 0)
dev.off()

png(filename = "./Figures/mature_celltype_Genotype_ridge2.png",width=30,height=15,units="cm",res=300)
RidgePlot(object = islets.integrated, cols=cluster_cols,#group.by ="Genotype", 
          features = mature[4:6], slot = 'counts', log = TRUE,
          pt.size = 0)
dev.off()

png(filename = "./Figures/mature_celltype_Genotype_ridge3.png",width=30,height=15,units="cm",res=300)
RidgePlot(object = islets.integrated, cols=cluster_cols,#group.by ="Genotype", 
          features = mature[7:9], slot = 'counts', log = TRUE,
          pt.size = 0)
dev.off()

png(filename = "./Figures/mature_celltype_Genotype_ridge4.png",width=30,height=15,units="cm",res=300)
RidgePlot(object = islets.integrated, cols=cluster_cols,#group.by ="Genotype", 
          features = mature[10:12], slot = 'counts', log = TRUE,
          pt.size = 0)
dev.off()


###### Cell Cycle markers #####
s.genes <- cc.genes$s.genes
s.genes <- paste("",s.genes,sep="")
g2m.genes <- cc.genes$g2m.genes
g2m.genes <- paste("",g2m.genes,sep="")

islets.integrated <- CellCycleScoring(islets.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)

# view cell cycle scores and phase assignments
head(islets.integrated[[]])

islets.integrated <- RunPCA(islets.integrated, features = c(s.genes, g2m.genes))
DimPlot(islets.integrated,group.by = "Phase")

saveRDS(islets.integrated, file = "./islets.integrated_final.rds")
#islets.integrated<-readRDS(file = "./islets.integrated_final.rds")

##### Run ALRA #####
#devtools::install_github('satijalab/seurat-wrappers')
library("SeuratWrappers")
islets.integrated
# Example 1: Simple usage, with automatic choice of k.
islets.integrated <- RunALRA(object = islets.integrated)

# visualize original and imputed values
islets.integrated <- NormalizeData(islets.integrated, assay = "alra")
features.plot <- endocrine[1]
DefaultAssay(islets.integrated) <- "RNA"
plot1 <- FeaturePlot(islets.integrated, features.plot, split.by="Genotype",ncol = 2)
DefaultAssay(islets.integrated) <- "alra"
plot2 <- FeaturePlot(islets.integrated, endocrine[1], ncol = 2, split.by="Genotype",cols = c("lightgrey", "red"))

png(filename = "./Figures/ins_ALRA.png",width=25,height=25,units="cm",res=300)
CombinePlots(plots=list(plot1,plot2),ncol=1)
dev.off()

##### FGSEA ####
library("org.Mm.eg.db")
library("fgsea")
library("gskb")

### Beta cells
idents1 <- c("0_KO","2_KO","3_KO")
idents2 <- c("0_WT","2_WT","3_WT")

b.cell.diffs_big <- FindMarkers(islets.integrated, ident.1 = idents1, 
                            ident.2 = idents2, 
                            assay="RNA",
                            slot = "data",
                            min.pct = 0.05, 
                            logfc.threshold = 0.05, verbose = T)

b.cell.diffs_big<-b.cell.diffs_big[order(b.cell.diffs_big$avg_logFC,b.cell.diffs_big$p_val_adj),]

ranks_beta <- b.cell.diffs_big$avg_logFC
names(ranks_beta) <- mapIdsList(x=org.Mm.eg.db, 
                                keys=row.names(b.cell.diffs_big),
                                keytype="SYMBOL", 
                                column="ENTREZID")

ranks_beta<-ranks_beta[!is.na(names(ranks_beta))]  

head(ranks_beta)

barplot(sort(ranks_beta, decreasing = T))
dev.off()

pathways_beta <- reactomePathways(names(ranks_beta))
fgseaRes_beta <- fgsea(pathways_beta, ranks_beta, maxSize=500)
head(fgseaRes_beta)

topPathwaysUp <- fgseaRes_beta[ES > 0][head(order(pval), n=20), pathway]
topPathwaysDown <- fgseaRes_beta[ES < 0][head(order(pval), n=20), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

png(filename = "./Figures/fgsea_reactome_beta.png",width=35,height=20,units="cm",res=300)
plotGseaTable(pathways_beta[topPathways], ranks_beta, fgseaRes_beta, gseaParam=0.5,
              colwidths = c(7, 1.5, 0.8, 0, 0.8),)
dev.off()

pathways_beta[topPathways]

#Activation of Matrix Metalloproteinases
mapIdsList(x=org.Mm.eg.db, 
           keys=c("16612","103964","22072","66473","22074"),
           keytype="ENTREZID", 
           column="SYMBOL")

#Digestion
mapIdsList(x=org.Mm.eg.db, 
           keys=c("18946","12613","109791","69060"),
           keytype="ENTREZID", 
           column="SYMBOL")

#Degradation of the extracellular matrix
mapIdsList(x=org.Mm.eg.db, 
           keys=c("12336","16612","103964","22072","66473","22074" ),
           keytype="ENTREZID", 
           column="SYMBOL")



#### Hallmark pathways
library("msigdbr")
h_gene_sets = msigdbr(species = "Mus musculus", category = "H")
head(h_gene_sets)
head(pathways_beta)

length(unique(h_gene_sets$gs_name))
pathways_hallmark<-as.list(unique(h_gene_sets$gs_name))
names(pathways_hallmark)<-pathways_hallmark
pathways_hallmark<-lapply(pathways_hallmark,function(l){
  as.character(h_gene_sets$entrez_gene[h_gene_sets$gs_name==l])
})

pathways_hallmark %>% 
  head() %>% 
  lapply(head)

fgseaRes.h.beta <- fgsea(pathways=pathways_hallmark, stats=ranks_beta)

fgseaResTidy.h.beta <- fgseaRes.h.beta %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy.h.beta %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

ggplot(fgseaResTidy.h.beta, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.1),show.legend = T) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

ggsave("./Figures/human_islets/hallmark_0.png",width=20,height=20,units="cm",dpi=300)

### acinar cells
idents1 <- c("12_KO","21_KO","20_KO","16_KO")
idents2 <- c("12_WT","21_WT","20_WT","16_WT")

acinar.diffs_big <-FindMarkers(islets.integrated, ident.1 = idents1, 
                               ident.2 = idents2, 
                               assay="RNA",
                               slot = "data",
                               min.pct = 0.05, 
                               logfc.threshold = 0.05, verbose = T)

acinar.diffs_big<-acinar.diffs_big[order(acinar.diffs_big$avg_logFC,acinar.diffs_big$p_val_adj),]

ranks_acinar <- acinar.diffs_big$avg_logFC
names(ranks_acinar) <- mapIdsList(x=org.Mm.eg.db, 
                                keys=row.names(acinar.diffs_big),
                                keytype="SYMBOL", 
                                column="ENTREZID")

ranks_acinar<-ranks_acinar[!is.na(names(ranks_acinar))]  

head(ranks_acinar)

barplot(sort(ranks_acinar, decreasing = T))
dev.off()

pathways_acinar <- reactomePathways(names(ranks_acinar))
fgseaRes_acinar <- fgsea(pathways_acinar, ranks_acinar, maxSize=500)
head(fgseaRes_acinar)

topPathwaysUp <- fgseaRes_acinar[ES > 0][head(order(pval), n=20), pathway]
topPathwaysDown <- fgseaRes_acinar[ES < 0][head(order(pval), n=20), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

png(filename = "./Figures/fgsea_reactome_acinar.png",width=35,height=20,units="cm",res=300)
plotGseaTable(pathways_acinar[topPathways], ranks_acinar, fgseaRes_acinar, gseaParam=0.5,
              colwidths = c(7, 1.5, 0.8, 0, 0.8),)
dev.off()

pathways_acinar[topPathways]

#Digestion
mapIdsList(x=org.Mm.eg.db, 
           keys=c("18946","12613","109791","69060"),
           keytype="ENTREZID", 
           column="SYMBOL")

#Degradation of the extracellular matrix
mapIdsList(x=org.Mm.eg.db, 
           keys=c("12336","16612","103964","22072","66473","22074" ),
           keytype="ENTREZID", 
           column="SYMBOL")


#### GO terms
pathways.GO <- gmtPathways("./human_islets/c5.all.v7.1.symbols.gmt")
pathways.GO %>% 
  head() %>% 
  lapply(head)

fgseaResGO.0 <- fgsea(pathways=pathways.GO, stats=ranks0, nperm=10000)

fgseaResTidy.GO.0 <- fgseaResGO.0 %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy.GO.0 %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

ggplot(fgseaResTidy.GO.0[fgseaResTidy.GO.0$pval<0.01,], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.1),show.legend = F) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GO terms NES from GSEA") + 
  theme_minimal()

ggsave("./Figures/human_islets/GO_0.png",width=20,height=20,units="cm",dpi=300)


ranks14 <- a.cell.diffs$avg_logFC[a.cell.diffs$cluster=="14"]
names(ranks14) <- a.cell.diffs$gene[a.cell.diffs$cluster=="14"]
head(ranks14)

#barplot(sort(ranks14, decreasing = T))


##### Perform GO Analysis using GOstats #####
library("biomartr")
library("GO.db")
library("org.Hs.eg.db")
library("GOstats")
library("categoryCompare")
library("RCy3")
library("pathview")
library("gage")
library("enrichplot")
library("clusterProfiler")
library("KEGG.db")

g2 <- mapIds(org.Hs.eg.db,
             keys=gsub("","",row.names(islets.integrated@assays$RNA@counts)),
             column="ENTREZID",
             keytype="SYMBOL",
             multiVals="first")
g2<-g2[!is.na(g2)]
head(g2)

x <- org.Hs.egGO
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])
got <- xx[[1]]
haveGo <- as.character(g2) %in% names(xx)
numNoGO <- sum(!haveGo) #616
g2<-g2[haveGo]

#Check individual duplicates by hand
universeIDs<-g2
dups<-data.frame(EntrezGene=universeIDs[duplicated(universeIDs)]) #no dups!
#dups$symbol<-mapIds(org.Hs.eg.db,
#                    keys=rownames(dups),
#                    column="SYMBOL",
#                    keytype="ENSEMBL",
#                    multiVals="first")
#dups

names(universeIDs)[duplicated(names(universeIDs))] #no duplicated ensembl ids

haveGo <- as.character(universeIDs) %in% names(xx)
numNoGO <- sum(!haveGo) #0
universeIDs[!haveGo] #sanity check
universeIDs<-universeIDs[haveGo]

##### Alpha cells #####
#Collect significantly different genes
a.cell.diffs$gene<-gsub("","",row.names(a.cell.diffs))
a.cell.diffs$ensembl<-mapIds(org.Hs.eg.db,
                             keys=a.cell.diffs$gene,
                             column="ENSEMBL",
                             keytype="SYMBOL",
                             multiVals="first")
a.cell.diffs$entrezid<-mapIds(org.Hs.eg.db,
                              keys=a.cell.diffs$gene,
                              column="ENTREZID",
                              keytype="SYMBOL",
                              multiVals="first")
a.cell.diffsGO<-a.cell.diffs[!is.na(a.cell.diffs$entrezid),]

geneLists<-list("Alpha-Like"=a.cell.diffsGO$entrezid)

geneLists<-lapply(geneLists,function(GL){
  GL<-list(genes=GL,universe=universeIDs,annotation="org.Hs.eg.db")
})

geneLists <- new("ccGeneList", geneLists, ccType=c("BP","CC","MF","KEGG")) #
geneLists

save.image("preenrich.RData")
enrichList <- ccEnrich(geneLists)

save(enrichList,file="enrichList.Rdata")
#load("enrichList.RData")
enrichList #check for lists with 0 nodes (p<0.05)

ego <- enrichGO(gene       = a.cell.diffsGO$entrezid[a.cell.diffsGO$p_val_adj<0.1],
                universe      = universeIDs,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.1,
                qvalueCutoff  = 0.1,
                readable      = TRUE)

enrichplot::dotplot(ego,showCategory = 30,orderBy="GeneRatio",color="qvalue")+
  ggtitle("Enriched Biological Processes in Alpha-like Cells from DV Grafts")
ggsave(filename = "./Figures/GO_BP_Cluster.png",width=25,height=15,units="cm")

egoMF <- enrichGO(gene       = a.cell.diffsGO$entrezid[a.cell.diffsGO$p_val_adj<0.1],
                  universe      = universeIDs,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
dotplot(ego,showCategory = 30,orderBy="GeneRatio",color="qvalue")+
  ggtitle("Enriched Molecular Functions in Alpha-like Cells from DV Grafts")

ggsave(filename = paste("./Figures/GO_MF_Cluster",clust,".png",sep=""),
       width=25,height=10,units="cm")



egoCC <- enrichGO(gene       = a.cell.diffsGO$entrezid[a.cell.diffsGO$p_val_adj<0.1],
                  universe      = universeIDs,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
dotplot(ego,showCategory = 30,orderBy="GeneRatio",color="qvalue")+
  ggtitle(paste("Enriched Cellular Components in Cluster ",clust,sep=""))
ggsave(filename = paste("./Figures/GO_CC_Cluster",clust,".png",sep=""),
       width=25,height=10,units="cm")


# kk <- enrichKEGG(gene         = gene,
#                  organism     = 'hsa',
#                  pvalueCutoff = 0.05)
save.image("postenrich.RData")

save(enrichList,file="enrichList.RData")
save(geneLists,file="geneLists.RData")
#save(GOLists,file="GOLists.RData")
#install.packages("VennDiagram")
#library("VennDiagram")
# ccOpts <- new("ccOptions", listNames=names(enrichList), 
#               compareNames=names(enrichList),
#               outType="None")
ccOptsBP <- new("ccOptions", listNames=names(enrichList$BP), 
                compareNames=names(enrichList$BP),
                outType="None")
ccOptsMF <- new("ccOptions", listNames=names(enrichList$MF), 
                compareNames=names(enrichList$MF),
                outType="None")
ccOptsCC <- new("ccOptions", listNames=names(enrichList$CC), 
                compareNames=names(enrichList$CC),
                outType="None")

# ccOptsCC <- new("ccOptions", listNames=c("Cluster 3","Cluster 4","Cluster 8"), 
#                 compareNames=c("Cluster 3","Cluster 4","Cluster 8"),
#                 #Don't compare between conditions - too many comparisons: could remove O, L, J, I
#                 outType="None")
# ccOptsMF <- new("ccOptions", listNames=c("Cluster 3","Cluster 4","Cluster 8"), 
#                 compareNames=c("Cluster 3","Cluster 4","Cluster 8"),
#                 #Don't compare between conditions - too many comparisons: could remove O, L, J, I
#                 outType="None")

ccOptsBP
save(ccOpts,file="ccOpts.RData")
save.image("preccCompare.RData")
#load("~/RProjects/CT_RNASeq/preccCompare.RData")
ccResultsBP <- ccCompare(ccEnrichResult=enrichList$BP,
                         ccOptions=ccOptsBP)

ccResultsMF <- ccCompare(ccEnrichResult=enrichList$MF,
                         ccOptions=ccOptsMF)

ccResultsCC <- ccCompare(ccEnrichResult=enrichList$CC,
                         ccOptions=ccOptsCC)

# ccResults <- ccCompare(ccEnrichResult=enrichList,
#                          ccOptions=ccOpts)

#save(ccResults,file="ccResults.RData")

load("./ccResults.RData")
load(file="./enrichListsub.RData")
#load("~/RProjects/ViaCyteRNASeq/gene.RData")

cwBP<-breakEdges(cwObject = ccResultsBP,0.8)
cwMF<-breakEdges(cwObject = ccResultsMF,0.2)
cwCC<-breakEdges(cwObject = ccResultsCC,0.2)

cw.BP <- ccOutCyt(cwBP,ccOptsBP)
cw.MF <- ccOutCyt(cwMF,ccOptsMF)
cw.CC <- ccOutCyt(cwCC,ccOptsCC)
#cw.KEGG<-ccOutCyt(ccResults$KEGG,ccOpts)


##### Gage Analysis with Pathview #####

#direction of fc, depends on levels(coldat$grp), the first level10
#taken as reference (or control) and the second one as experiment.
kg.hsa <-kegg.gsets(species = "hsa", id.type = "kegg", check.new=FALSE)
kegg.hsa.gs<-kg.hsa$kg.sets[kg.hsa$sigmet.idx]

fc<-a.cell.diffsGO$avg_logFC
names(fc)<-as.character(a.cell.diffsGO$entrezid)
fc.kegg.p <- gage(fc, gsets = kegg.hsa.gs, ref = NULL, samp = NULL)
sel <- fc.kegg.p$greater[, "q.val"] < 0.1 &
  !is.na(fc.kegg.p$greater[, "q.val"])
path.ids <- rownames(fc.kegg.p$greater)[sel]
sel.l <- fc.kegg.p$less[, "q.val"] < 0.1 &
  !is.na(fc.kegg.p$less[,"q.val"])
path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)
pv.out.list <- sapply(path.ids2, function(pid) pathview(gene.data =  ct.fc, pathway.id = pid,
                                                        low = list(gene = "#4D9221", cpd = "blue"), 
                                                        mid = list(gene = "#F6F6F5", cpd= "gray"), 
                                                        high = list(gene = "#C72381", cpd = "yellow"),
                                                        species = "dre", out.suffix=paste("Cluster",clust,sep=""),
                                                        kegg.native = T,cex=0.25,
                                                        kegg.dir = paste("./Figures/KEGG/Cluster",clust,sep="")))

download.kegg("04010","dre",kegg.dir="./Figures/KEGG/Cluster0")

##### All clusters #####
#Collect significantly different genes
rat_markers$ensembl<-mapIds(org.Hs.eg.db,
                            keys=rat_markers$gene,
                            column="ENSEMBL",
                            keytype="SYMBOL",
                            multiVals="first")
rat_markers$entrezid<-mapIds(org.Hs.eg.db,
                             keys=rat_markers$gene,
                             column="ENTREZID",
                             keytype="SYMBOL",
                             multiVals="first")
rat_markersGO<-rat_markers[!is.na(rat_markers$entrezid),]

geneLists<-split(rat_markersGO$entrezid,rat_markersGO$cluster)
#geneLists<-geneLists[!names(geneLists)%in%c("5","6")] #Clusters 5 and 6 have no enriched processes

geneLists<-lapply(geneLists,function(GL){
  GL<-list(genes=GL,universe=universeIDs,annotation="org.Hs.eg.db")
})

save.image("preenrich.RData")
enrichList <- ccEnrich(geneLists)

save(enrichList,file="enrichList.Rdata")
#load("enrichList.RData")
enrichList #check for lists with 0 nodes (p<0.05)

save(enrichList,file="enrichList.Rdata")
#load("enrichList.RData")
enrichList #check for lists with 0 nodes (p<0.05)
# ccZerosBP<-c("Cluster 1","Cluster 7")
# ccZerosCC<-c("Cluster 2","Cluster 7")
# ccZerosMF<-c("Cluster 7")
#ccZerosKEGG<-c("Cluster 6")

#pvalueCutoff(enrichList$BP) <- 0.001 #If too many nodes
# enrichList$BP<-enrichList$BP[!(names(enrichList$BP)%in%ccZerosBP)]
# enrichList$CC<-enrichList$CC[!(names(enrichList$CC)%in%ccZerosCC)]
# enrichList$MF<-enrichList$MF[!(names(enrichList$MF)%in%ccZerosMF)]
# enrichList$KEGG<-enrichList$KEGG[!(names(enrichList$KEGG)%in%ccZerosKEGG)]

# BPList<-lapply(enrichList$BP,summary)
# MFList<-lapply(enrichList$MF,summary)
# CCList<-lapply(enrichList$CC,summary)
# KEGGList<-lapply(enrichList$KEGG,summary)
# GOLists<-list(BPList=BPList,
#               MFList=MFList,
#               CCList=CCList)

#clustBP<-names(enrichList$BP[!(names(enrichList$BP)%in%ccZerosBP)])
clustBP<-names(enrichList$BP)
clustBP<-gsub("Cluster ","",clustBP)

#clustCC<-names(enrichList$CC[!(names(enrichList$CC)%in%ccZerosCC)])
clustCC<-names(enrichList$CC)
clustCC<-gsub("Cluster ","",clustCC)

#clustMF<-names(enrichList$MF[!(names(enrichList$MF)%in%ccZerosMF)])
clustMF<-names(enrichList$MF)
clustMF<-gsub("Cluster ","",clustMF)

clustKEGG<-names(enrichList$KEGG)
#clustKEGG<-names(enrichList$KEGG)
clustKEGG<-gsub("Cluster ","",clustKEGG)

for (clust in clustBP ){
  ego <- enrichGO(gene       = rat_markersGO$entrezid[rat_markersGO$cluster==clust&
                                                        rat_markersGO$p_val_adj<0.05],
                  universe      = universeIDs,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  dotplot(ego,showCategory = 30,orderBy="GeneRatio",color="qvalue")+
    ggtitle(paste("Enriched Biological Processes in Cluster ",clust,sep=""))
  ggsave(filename = paste("./Figures/GO/BP/GO_BP_Cluster",clust,".png",sep=""),
         width=25,height=15,units="cm")
}

for (clust in clustMF ){
  ego <- enrichGO(gene       = rat_markersGO$entrezid[rat_markersGO$cluster==clust&
                                                        rat_markersGO$p_val_adj<0.05],
                  universe      = universeIDs,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  dotplot(ego,showCategory = 30,orderBy="GeneRatio",color="qvalue")+
    ggtitle(paste("Enriched Molecular Functions in Cluster ",clust,sep=""))
  ggsave(filename = paste("./Figures/GO/MF/GO_MF_Cluster",clust,".png",sep=""),
         width=25,height=10,units="cm")
}

for (clust in clustCC ){
  ego <- enrichGO(gene       = rat_markersGO$entrezid[rat_markersGO$cluster==clust&
                                                        rat_markersGO$p_val_adj<0.05],
                  universe      = universeIDs,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  dotplot(ego,showCategory = 30,orderBy="GeneRatio",color="qvalue")+
    ggtitle(paste("Enriched Cellular Components in Cluster ",clust,sep=""))
  ggsave(filename = paste("./Figures/GO/CC/GO_CC_Cluster",clust,".png",sep=""),
         width=25,height=10,units="cm")
}
gc()

# kk <- enrichKEGG(gene         = gene,
#                  organism     = 'hsa',
#                  pvalueCutoff = 0.05)
save.image("postenrich.RData")

#save(GOLists,file="GOLists.RData")
#install.packages("VennDiagram")
#library("VennDiagram")
# ccOpts <- new("ccOptions", listNames=names(enrichList), 
#               compareNames=names(enrichList),
#               outType="None")
ccOptsBP <- new("ccOptions", listNames=names(enrichList$BP), 
                compareNames=names(enrichList$BP),
                outType="None")
ccOptsMF <- new("ccOptions", listNames=names(enrichList$MF), 
                compareNames=names(enrichList$MF),
                outType="None")
ccOptsCC <- new("ccOptions", listNames=names(enrichList$CC), 
                compareNames=names(enrichList$CC),
                outType="None")

# ccOptsCC <- new("ccOptions", listNames=c("Cluster 3","Cluster 4","Cluster 8"), 
#                 compareNames=c("Cluster 3","Cluster 4","Cluster 8"),
#                 #Don't compare between conditions - too many comparisons: could remove O, L, J, I
#                 outType="None")
# ccOptsMF <- new("ccOptions", listNames=c("Cluster 3","Cluster 4","Cluster 8"), 
#                 compareNames=c("Cluster 3","Cluster 4","Cluster 8"),
#                 #Don't compare between conditions - too many comparisons: could remove O, L, J, I
#                 outType="None")

ccOptsBP
save(ccOpts,file="ccOpts.RData")
save.image("preccCompare.RData")
#load("~/RProjects/CT_RNASeq/preccCompare.RData")
ccResultsBP <- ccCompare(ccEnrichResult=enrichList$BP,
                         ccOptions=ccOptsBP)

ccResultsMF <- ccCompare(ccEnrichResult=enrichList$MF,
                         ccOptions=ccOptsMF)

ccResultsCC <- ccCompare(ccEnrichResult=enrichList$CC,
                         ccOptions=ccOptsCC)

# ccResults <- ccCompare(ccEnrichResult=enrichList,
#                          ccOptions=ccOpts)

#save(ccResults,file="ccResults.RData")

load("./ccResults.RData")
load(file="./enrichListsub.RData")
#load("~/RProjects/ViaCyteRNASeq/gene.RData")

cwBP<-breakEdges(cwObject = ccResultsBP,0.8)
cwMF<-breakEdges(cwObject = ccResultsMF,0.2)
cwCC<-breakEdges(cwObject = ccResultsCC,0.2)

cw.BP <- ccOutCyt(cwBP,ccOptsBP)
cw.MF <- ccOutCyt(cwMF,ccOptsMF)
cw.CC <- ccOutCyt(cwCC,ccOptsCC)
#cw.KEGG<-ccOutCyt(ccResults$KEGG,ccOpts)


##### Gage Analysis with Pathview #####

#direction of fc, depends on levels(coldat$grp), the first level10
#taken as reference (or control) and the second one as experiment.
kg.hsa <-kegg.gsets(species = "hsa", id.type = "kegg", check.new=FALSE)
kegg.hsa.gs<-kg.hsa$kg.sets[kg.hsa$sigmet.idx]

fc<-a.cell.diffsGO$avg_logFC
names(fc)<-as.character(a.cell.diffsGO$entrezid)
fc.kegg.p <- gage(fc, gsets = kegg.hsa.gs, ref = NULL, samp = NULL)
sel <- fc.kegg.p$greater[, "q.val"] < 0.1 &
  !is.na(fc.kegg.p$greater[, "q.val"])
path.ids <- rownames(fc.kegg.p$greater)[sel]
sel.l <- fc.kegg.p$less[, "q.val"] < 0.1 &
  !is.na(fc.kegg.p$less[,"q.val"])
path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)
pv.out.list <- sapply(path.ids2, function(pid) pathview(gene.data =  ct.fc, pathway.id = pid,
                                                        low = list(gene = "#4D9221", cpd = "blue"), 
                                                        mid = list(gene = "#F6F6F5", cpd= "gray"), 
                                                        high = list(gene = "#C72381", cpd = "yellow"),
                                                        species = "dre", out.suffix=paste("Cluster",clust,sep=""),
                                                        kegg.native = T,cex=0.25,
                                                        kegg.dir = paste("./Figures/KEGG/Cluster",clust,sep="")))

##### Look at differences in insulin positive cells between Genotypes #####
Idents(islets.integrated) <- "Genotype"
ins.pos<-WhichCells(islets.integrated,expression = `INS`>3250,slot="counts" )
#or by cluster, ideally

png(filename = "./Figures/heatmap_inspos_goi.png",width=15,height=20,units="cm",res=300)
DoHeatmap(object = islets.integrated, cells=ins.pos, assay="RNA",group.colors = rev(Genotype_cols), 
          features = c(endocrine,biosynthesis,ECL,endothelial,
                       immune,mature,metab,stromal), angle=0)
dev.off()

betagenes<-as.character(row.names(rat_markers)[rat_markers$cluster=="β-like"])

png(filename = "./Figures/rat_heatmap_inspos.png",width=15,height=10,units="cm",res=300)
DoHeatmap(object = islets.integrated, cells=ins.pos,group.colors = rev(Genotype_cols), 
          features = top10$gene[top10$cluster=="β-like"], angle=0)
dev.off()

islets.integrated$celltype[ip_ins.pos]
