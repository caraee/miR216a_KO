##### Set-up #####
library("here")
library("Biobase")
library("dplyr")
library("tidyr")
library("ggplot2")
library("RColorBrewer")
library("extrafont")
library("GOstats")
library("DESeq2")
library("Rgraphviz")
library("GO.db")
library("KEGG.db")
library("KEGGREST")
library("genefilter")
library("org.Mm.eg.db")
library("igraph")
library("categoryCompare")
library("RCy3")
library("pathview")
library("gage")
#devtools::install_github("gsimchoni/kandinsky")
save.image("~/RProjects/Suheda_RNASeq/pre-GO_workspace.RData")

##### Filter genes for GO #####
#Keep only most variable genes
iqrCutoff <- 0.5
fc_filtered_IQR <- apply(counts(ddsTxi), 1, IQR)
selected <- fc_filtered_IQR > iqrCutoff
ddsTxiGO <- ddsTxi[selected, ]
g2 <- mapIds(org.Mm.eg.db,
             keys=row.names(ddsTxiGO),
             column="ENTREZID",
             keytype="ENSEMBL",
             multiVals="first")
g2<-g2[!is.na(g2)]
head(g2)
ddsTxiGO<-ddsTxiGO[names(g2),]
agg<-aggregate(counts(ddsTxiGO), by=list(g2), FUN=mean)
row.names(agg)<-as.character(agg$Group.1)
agg$Group.1<-NULL

#Remove genes with no GO mapping
x <- org.Mm.egGO
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])
got <- xx[[1]]
haveGo <- as.character(g2) %in% names(xx)
numNoGO <- sum(!haveGo) #235
ddsTxiGO <- ddsTxiGO[haveGo, ]
g2<-g2[haveGo]

#Check individual duplicates by hand
universeIDs<-row.names(agg)
dups<-data.frame(EntrezGene=universeIDs[duplicated(universeIDs)]) #no dups!
#dups$symbol<-mapIds(org.Hs.eg.db,
#                    keys=rownames(dups),
#                    column="SYMBOL",
#                    keytype="ENSEMBL",
#                    multiVals="first")
#dups

names(universeIDs)[duplicated(names(universeIDs))] #no duplicated ensembl ids

haveGo <- as.character(universeIDs) %in% names(xx)
numNoGO <- sum(!haveGo) #234
universeIDs[!haveGo] #sanity check
universeIDs<-universeIDs[haveGo]

##### Perform GO Analysis using GOstats #####
library("biomartr")

#Collect significantly different genes
universeIDs<-row.names(agg)
#results<-as.data.frame(results)
#names(results)<-c(names(resListsub[[1]]),"ensembl_gene_id","hgnc_symbol","gene_biotype","entrezgene")
genes2<-as.character(unique(results$entrezgene[which(results$padj<0.05)]))

names(genes2)<-mapIds(org.Mm.eg.db,
                       keys=genes2,
                       column="ENSEMBL",
                       keytype="ENTREZID",
                       multiVals="first")
genes2<-genes2[!is.na(genes2)]
genes2<-genes2[genes2%in%universeIDs]
genes<-list(genes=genes2,universe=universeIDs,annotation="org.Mm.eg.db")
genes<-list(WTvsKO=genes)

results_GO<-results[results$entrezgene%in%genes2,]
write.csv(results_GO,file="results_GO.csv",row.names=T, quote=F)

geneList <- new("ccGeneList", genes, ccType=c("BP","CC","MF","KEGG"))
geneList
enrichList <- ccEnrich(geneList)
save(enrichList,file="enrichList.Rdata")
#load("~/RProjects/Suheda_RNASeq/enrichList.RData")
enrichList #check for lists with 0 nodes (p<0.05)
#ccNonZero<-c("Cond1vs3","Cond1vs6","Cond1vsHI","Cond2vs3","Cond2vs6","Cond2vsHI",
#             "Cond3vs5","Cond3vs6","Cond3vsHI","Cond5vs6","Cond5vsHI","Cond6vsHI")
#pvalueCutoff(enrichList$BP) <- 0.001 #If too many nodes

BPList<-lapply(enrichList$BP,summary)
MFList<-lapply(enrichList$MF,summary)
CCList<-lapply(enrichList$CC,summary)
KEGGList<-lapply(enrichList$KEGG,summary)
GOLists<-list(BPList=BPList,
              MFList=MFList,
              CCList=CCList,
              KEGGList=KEGGList)

xx <- as.list(org.Mm.egGO2ALLEGS)

#rownames(vsd.fast)<-filteredIDs
GOLists2<-lapply(names(GOLists),function(go){
  go<-lapply(names(GOLists[[go]]),function(nm){
    df<-GOLists[[go]][[nm]]
    df<-df[df$FDR<0.05,]
    if (nrow(df[df$Size<=500,])>5) n=5 else n=nrow(df[df$Size<=500,])
    if (n>0){
      df<-df[df$Size<=500,]%>%.[order(.$OddsRatio,decreasing = T),]%>%.[1:n,]
      df$Label<-AnnotationDbi::select(GO.db,keys=df$ID,
                                      columns="TERM",keytype="GOID")
      GOgenes<-sapply(df$ID,function(goid){
        genes<-unlist(xx[goid])
        genes<-genes[genes%in%geneList[[nm]]$genes]
        genes<-mapIds(org.Mm.eg.db,
                      keys=genes,column="ENSEMBL",
                      keytype="ENTREZID",multiVals="first")})
      lapply(names(GOgenes),function(goid){
        goi<-GOgenes[[goid]]
        goi<-goi[goi%in%rownames(ddsTxiGO)]
        padj<-results[goi,]$padj
        names(padj)<-goi
        padj<-padj[unique(names(padj))]
        padj<-padj[order(padj)]
        if (length(padj)>12) padj<-padj[1:12]
        label<-df$Label$TERM[df$ID==goid]
        tcounts<-t(counts(ddsTxiGO,normalized=T,replaced=T)[rownames(ddsTxiGO)%in%names(padj),]) %>%
          merge(colData(ddsTxiGO), ., by="row.names") %>%
          gather(gene, expression, (ncol(.)-length(names(padj))+1):ncol(.))
        tcounts$symbol<-mapIds(org.Mm.eg.db,
                               keys=names(padj),column="SYMBOL",
                               keytype="ENSEMBL",multiVals="first")
        main<-"WT vs KO"
        cond<-unlist(strsplit(nm,""))[5]
        p<- ggplot(tcounts, aes(Genotype, expression, fill=Genotype))
        p+geom_boxplot() + 
          facet_wrap(~symbol, scales="free_y") + 
          #scale_y_continuous(trans="log2")+
          scale_color_manual(name="Genotype",values = ann_colors$Genotype)+
          scale_fill_manual(name="Genotype",values = ann_colors$Genotype)+
          #guides(colour=guide_legend(ncol=2),fill=guide_legend(ncol=2))+
          xlab("")+
          labs(x="",y="Expression\n(Normalized counts, outliers replaced)", 
               fill="Genotype", title=paste(main,":\n",label,sep=""))+
          theme_bw()+
          theme(plot.title = element_text(family = "Open Sans Light", color="#666666",  size=20, hjust=0)) +
          theme(axis.title = element_text(family = "Roboto Light", color="#666666", size=18))+
          theme(axis.text.x = element_blank())+
          theme(legend.text = element_text(family = "Roboto Light")) #+
        #theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))
        ggsave(paste(go,nm,goid,".png",sep=""),path="./Figures/GO",
               width = 20, height = 15, units = "cm")
        })
    }
  })
})

df<-GOLists[["BPList"]][["WTvsKO"]]
write.csv(df,file="BP_treemap.csv")

df<-GOLists[["CCList"]][["WTvsKO"]]
write.csv(df,file="CC_treemap.csv")

df<-GOLists[["MFList"]][["WTvsKO"]]
write.csv(df,file="MF_treemap.csv")

save(enrichList,file="enrichList.RData")
#save(geneLists,file="geneLists.RData")
save(GOLists,file="GOLists.RData")
#install.packages("VennDiagram")
library("VennDiagram")
ccOpts <- new("ccOptions", listNames="WTvsKO", 
              compareNames=names(geneList), #Don't compare between conditions - too many comparisons: could remove O, L, J, I
              outType=c("text","html","rcy3"))
ccOpts
save(ccOpts,file="ccOpts.RData")
ccResults <- ccCompare(enrichList,ccOpts)

save(ccResults,file="ccResults.RData")

cwBP<-breakEdges(cwObject = ccResults$BP,0.6)
cwMF<-breakEdges(cwObject = ccResults$MF,0.2)
cwCC<-breakEdges(cwObject = ccResults$CC,0.2)

cw.BP <- ccOutCyt(cwBP,ccOpts)
cw.MF <- ccOutCyt(cwMF,ccOpts)
cw.CC <- ccOutCyt(cwCC,ccOpts)
cw.KEGG<-ccOutCyt(ccResults$KEGG,ccOpts)




##### Gage Analysis with Pathview #####
cnts.kegg.p <- gage(cnts.norm, gsets = kegg.gs, ref = ref.idx,
                    samp = samp.idx, compare ="unpaired")
cnts.d= cnts.norm[, samp.idx]-rowMeans(cnts.norm[, ref.idx])
sel <- cnts.kegg.p$greater[, "q.val"] < 0.1 &
  !is.na(cnts.kegg.p$greater[,"q.val"])
path.ids <- rownames(cnts.kegg.p$greater)[sel]
sel.l <- cnts.kegg.p$less[, "q.val"] < 0.1 &!is.na(cnts.kegg.p$less[,"q.val"])
path.ids.l <- rownames(cnts.kegg.p$less)[sel.l]
path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)
pv.out.list <- sapply(path.ids2, function(pid) pathview(gene.data = cnts.d, pathway.id = pid,
                                                        species = "hsa"))

grp.idx <- rep(c("knockdown", "control"), each=4)
coldat=DataFrame(grp=factor(grp.idx))
dds <- DESeqDataSetFromMatrix(cnts, colData=coldat, design = ~ grp)
dds <- DESeq(dds)
deseq2.res <- results(dds)
#direction of fc, depends on levels(coldat$grp), the first level10
  #taken as reference (or control) and the second one as experiment.

resWTvsKO.fc<-resWTvsKO$log2FoldChange
names(resWTvsKO.fc)<-as.character(resWTvsKO$entrezgene)
exp.fc<-resWTvsKO.fc
out.suffix<-"deseq2"

kg.mmu<-kegg.gsets(species = "mmu", id.type = "kegg", check.new=FALSE)
kegg.mmu.gs<-kg.mmu$kg.sets[kg.mmu$sigmet.idx]

fc.kegg.p <- gage(exp.fc, gsets = kegg.mmu.gs, ref = NULL, samp = NULL)

sel <- fc.kegg.p$greater[, "q.val"] < 0.1 &
  !is.na(fc.kegg.p$greater[, "q.val"])
path.ids <- rownames(fc.kegg.p$greater)[sel]
sel.l <- fc.kegg.p$less[, "q.val"] < 0.1 &
  !is.na(fc.kegg.p$less[,"q.val"])
path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)
pv.out.list <- sapply(path.ids2, function(pid) pathview(gene.data =  exp.fc, pathway.id = pid,
                                                        low = list(gene = "#4D9221", cpd = "blue"), 
                                                        mid = list(gene = "#F6F6F5", cpd= "gray"), 
                                                        high = list(gene = "#C72381", cpd = "yellow"),
                                                             species = "mmu", out.suffix=out.suffix,
                                                        kegg.native = T,cex=0.25,
                                                        kegg.dir = "./Figures/KEGG"))

