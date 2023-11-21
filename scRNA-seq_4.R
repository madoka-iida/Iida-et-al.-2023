#The scRNA-seq dataset was processed and analyzed with Seurat version 3.1.460.
# Load libraries
library(Seurat) #version 3.1.460
library(cowplot)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(rvcheck)
library(org.Mm.eg.db)
BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
library(gridExtra)
library(grid)

#Loading the dataset
W3.data <- Read10X_h5("W3 filtered_feature_bc_matrix-2.h5")
A3.data <- Read10X_h5("A3 filtered_feature_bc_matrix-2.h5")

#Creating the seurat object
W3 <- CreateSeuratObject(counts = W3.data,project = "W3", min.cell = 3)
A3 <- CreateSeuratObject(counts = A3.data,project = "A3", min.cell = 3)

#Store in a group named WA3
W3$WA3 <- "W3"
A3$WA3 <- "A3"

#QC and selecting cells for further analysis
W3[["percent.mt"]] <- PercentageFeatureSet(W3, pattern = "mt-")
W3 <- subset(W3, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 1)
A3[["percent.mt"]] <- PercentageFeatureSet(A3, pattern = "mt-")
A3 <- subset(A3, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 1)

#Make a list
WA3.list <- NULL
WA3.list <- list(W3, A3)
names(WA3.list) <- c("W3", "A3")

#Normalize and identify variable features for each dataset independently
WA3.list <- lapply(X = WA3.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)  
})

#Select features that are repeatedly variable across datasets for integration run PCA on each
features <- SelectIntegrationFeatures(object.list = WA3.list)
WA3.list <- lapply(X = WA3.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

#Perform integration
WA3.anchors <- FindIntegrationAnchors(object.list = WA3.list, anchor.features = features, reduction = "rpca", k.anchor =2 )
WA3.combined <- IntegrateData(anchorset = WA3.anchors)
DefaultAssay(WA3.combined) <- "integrated"

#Visualization and clustering (tSNE)
WA3.combined <- ScaleData(WA3.combined, verbose = FALSE)
WA3.combined <- RunPCA(WA3.combined, npcs = 25, verbose = FALSE)
WA3.combined <- RunTSNE(WA3.combined, reduction = "pca", dims = 1:25)
WA3.combined  <- FindNeighbors(WA3.combined, graph.name = "WA3test", reduction = "pca", dims = 1:25)
WA3.combined  <- FindClusters(WA3.combined, graph.name = "WA3test", resolution = 0.8)

#Visualization and clustering (UMAP)
WA3.combined <- ScaleData(WA3.combined, verbose = FALSE)
WA3.combined <- RunPCA(WA3.combined, npcs = 25, verbose = FALSE)
WA3.combined <- RunUMAP(WA3.combined, reduction = "pca", dims = 1:25)
WA3.combined  <- FindNeighbors(WA3.combined, graph.name = "WA3test", reduction = "pca", dims = 1:25)
WA3.combined  <- FindClusters(WA3.combined, graph.name = "WA3test", resolution = 0.8)

DefaultAssay(WA3.combined) <- "RNA" 

##DotPlot 
features <- c("Pdgfra", "C1ql1", "Tnr", "Ust", "Cnksr3", "Ptgds", "Mobp", "Mog", "Mag", "Mbp", "Plp1", "Celf4", "Rbfox3", "Slc17a6", "Slc6a1", "Slc6a5", "Gad2", "Chat","Slc1a2", "Aqp4", "Ly86", "Runx1", "Vtn", "Abcc9", "Flt1","Dnah12")
DotPlot(WA3.combined, features = features) 

#Exclude clusters containing doublets or background
WA3.combined.subset <- subset(WA3.combined,  idents = c("11","19", "21"), invert=TRUE)
saveRDS(WA3.combined.subset,"~/Downloads/R/230908 3w RPCA_before_annotation.rds") 

#Heatmap
WA3.combined.marker <- FindAllMarkers(WA3.combined.subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
WA3.combined.marker %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top5 <- WA3.combined.marker %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(WA3.combined.subset, features = top5$gene) + NoLegend()

#DotPlot  
features <- c("Pdgfra", "C1ql1", "Tnr", "Ust", "Cnksr3", "Ptgds", "Mobp", "Mog", "Mag", "Mbp", "Plp1", "Celf4", "Rbfox3", "Slc17a6", "Slc6a1", "Slc6a5", "Gad2", "Chat","Slc1a2", "Aqp4", "Ly86", "Runx1", "Vtn", "Abcc9", "Flt1","Dnah12")
DotPlot(WA3.combined.subset, features = features) 


#Label clusters 
new.cluster.ids <- c("Oligodendrocytes","Oligodendrocytes","Inhibitory neurons","Oligodendrocytes",  
                     "Astrocytes", "Inhibitory neurons","Oligodendrocytes","OPC",
                     "Inhibitory neurons","Excitatory neurons","OL progenitors", "Excitatory neurons",
                     "Microglia","Excitatory neurons","Excitatory neurons","Inhibitory neurons",
                     "Excitatory neurons","Excitatory neurons","Pericytes", "Inhibitory neurons")
names(new.cluster.ids) <- levels(WA3.combined.subset)
WA3.combined.subset <- RenameIdents(WA3.combined.subset, new.cluster.ids)
DimPlot(WA3.combined.subset , reduction = "tsne", label = TRUE, split.by="WA3")
DimPlot(WA3.combined.subset, reduction = "tsne", label = FALSE, 
        cols = c("Oligodendrocytes"="#F8766D","Inhibitory neurons" ="#CD9600", "Astrocytes"="#7CAE00", "Excitatory neurons"="#0CB702", "OPC"="#ED68ED", "OL progenitors"="#00A9FF", "Microglia"="#00BFC4", "Pericytes"="#8494FF" ), 
        order = c("Pericytes","Microglia","OL progenitors","OPC","Astrocytes","Excitatory neurons", "Inhibitory neurons","Oligodendrocyte"))  ###採用Fig.2B
DimPlot(WA3.combined.subset, reduction = "umap", label = FALSE,
        cols = c("Oligodendrocytes"="#F8766D","Inhibitory neurons" ="#CD9600", "Astrocytes"="#7CAE00", "Excitatory neurons"="#0CB702", "OPC"="#ED68ED", "OL progenitors"="#00A9FF", "Microglia"="#00BFC4", "Pericytes"="#8494FF" ),
        order = c("Pericytes","Microglia","OL progenitors","OPC","Astrocytes","Excitatory neurons", "Inhibitory neurons","Oligodendrocyte"))  ###採用Fig.2B
table(Idents(WA3.combined.subset), WA3.combined.subset$WA3) 


#FeaturePlot
oligo <- subset(x = WA3.combined.subset, idents=c("Oligodendrocytes"))
features <- c("Asic2")  
FeaturePlot(object = oligo, features = features,split.by="WA3", reduction = "tsne",pt.size = 0.2,   min.cutoff = 2, max.cutoff = 3.5)
features <- c("Fam155a")
FeaturePlot(object = oligo, features = features,split.by="WA3",reduction = "tsne", pt.size = 0.2,   min.cutoff = 1.5, max.cutoff =3.5)


WA3.combined.subset <- readRDS("~/Downloads/R/230908 3w RPCA_before_annotation.rds")
#DEGs of OLs clusters
WA3.combined.subset$cluster_sample <- paste(WA3.combined.subset$seurat_clusters, WA3.combined.subset$orig.ident, sep = "_")
Idents(WA3.combined.subset) <- "cluster_sample"
DEGoligo_WA3 <- FindMarkers(WA3.combined.subset, ident.1 = c("0_A3","1_A3","3_A3","6_A3"), ident.2 = c("0_W3","1_W3","3_W3","6_W3"), logfc.threshold = 0.1)

#GO enrichment analysis of 18 genes elevated in AR-97Q mice
DEGoligo_WA3_0.05_0.2 <- DEGoligo_WA3[DEGoligo_WA3$p_val_adj < 0.05
                                      & DEGoligo_WA3$avg_log2FC>
                                        0.2,]
DEG.gene_SYMBOLs <- rownames(DEGoligo_WA3_0.05_0.2)
Mm<- org.Mm.eg.db
DEG.gene_IDs <- AnnotationDbi::select(Mm,
                                      keys=DEG.gene_SYMBOLs,
                                      columns = c(
                                        "SYMBOL",
                                        "ENTREZID"
                                      ),
                                      keytype= "SYMBOL"  )$ENTREZID

#MF
ego <- enrichGO(gene = DEG.gene_IDs,
                ont ="MF",
                OrgDb='org.Mm.eg.db')
dotplot(ego, showCategory=5)  

#BP
ego <- enrichGO(gene = DEG.gene_IDs,
                ont ="BP",
                OrgDb='org.Mm.eg.db')
dotplot(ego, showCategory=5) 


#GO enrichment analysis of 67 genes downregulated in AR-97Q mice
DEGoligo_WA3_0.05_0.2 <- DEGoligo_WA3[DEGoligo_WA3$p_val_adj < 0.05
                                      & DEGoligo_WA3$avg_log2FC< -0.2
                                      ,]

DEG.gene_SYMBOLs <- rownames(DEGoligo_WA3_0.05_0.2)
Mm<- org.Mm.eg.db
DEG.gene_IDs <- AnnotationDbi::select(Mm,
                                      keys=DEG.gene_SYMBOLs,
                                      columns = c(
                                        "SYMBOL",
                                        "ENTREZID"
                                      ),
                                      keytype= "SYMBOL"  )$ENTREZID
#MF
ego <- enrichGO(gene = DEG.gene_IDs,
                ont ="MF",
                OrgDb='org.Mm.eg.db')
dotplot(ego, showCategory=5) 
#BP
ego <- enrichGO(gene = DEG.gene_IDs,
                ont ="BP",
                OrgDb='org.Mm.eg.db')
dotplot(ego, showCategory=5)  

#Volcano plot
WA3.combined.subset <- readRDS("~/Downloads/R/230908 3w RPCA_before_annotation.rds")
DEGoligo_WA3 <- FindMarkers(WA3.combined.subset, ident.1 = c("0_A3","1_A3","3_A3","6_A3"), ident.2 = c("0_W3","1_W3","3_W3","6_W3"), logfc.threshold = 0.1)

#Annotate the Ensembl gene IDs to gene symbols:
ens <- rownames(DEGoligo_WA3)
symbols <- mapIds(org.Mm.eg.db, keys = ens,
                  column = c('SYMBOL'), keytype = 'SYMBOL')
symbols <- symbols[!is.na(symbols)]
symbols <- symbols[match(rownames(DEGoligo_WA3), names(symbols))]
rownames(DEGoligo_WA3) <- symbols
keep <- !is.na(rownames(DEGoligo_WA3))
DEGoligo_WA3  <- DEGoligo_WA3 [keep,]

EnhancedVolcano(DEGoligo_WA3 , 
                lab = rownames(DEGoligo_WA3),
                x ="avg_log2FC", 
                y ="p_val_adj",
                pCutoff = 10e-2, 
                FCcutoff = 0.2,
                pointSize = 1.0,selectLab = c('Asic2', "Rbfox1", "Kcnip4", "Fam155a","Camk1d","Gm47283", "Il31ra", "Fmn1", "Samd4", "Bach2" ),
                labSize = 5.0,
                drawConnectors=TRUE, widthConnectors=0.5) 

#Reanalyze the OL cluster
#Visualization and clustering (tSNE)
all.genes <- rownames(WA3.combined.OL)
WA3.combined.OL <- ScaleData(WA3.combined.OL, features = all.genes)
WA3.combined.OL <- FindVariableFeatures(object = WA3.combined.OL)
WA3.combined.OL <- RunPCA(WA3.combined.OL, features = VariableFeatures(object = WA3.combined.OL))
WA3.combined.OL <- RunTSNE(WA3.combined.OL , reduction = "pca", dims = 1:25)
WA3.combined.OL  <- FindNeighbors(WA3.combined.OL, graph.name = "WA3test", reduction = "pca", dims = 1:25)
WA3.combined.OL  <- FindClusters(WA3.combined.OL, graph.name = "WA3test", resolution = 0.8)
p1 <- DimPlot(WA3.combined.OL , reduction = "tsne", group.by = "WA3")
plot_grid(p1)　

#Visualization and clustering (UMAP)
WA3.combined.OL <- RunUMAP(WA3.combined.OL , reduction = "pca", dims = 1:25)
WA3.combined.OL  <- FindNeighbors(WA3.combined.OL, graph.name = "WA3test", reduction = "pca", dims = 1:25)
WA3.combined.OL  <- FindClusters(WA3.combined.OL, graph.name = "WA3test", resolution = 0.8)
p1 <- DimPlot(WA3.combined.OL , reduction = "umap", group.by = "WA3")
plot_grid(p1)　

#Subset the oligodendrocyte cluster　
dataoligo <- subset(x = WA3.combined.OL, idents = c("0", "1","2","3","4","5"))
table(Idents(dataoligo), dataoligo$WA3) 
DimPlot(dataoligo, label = FALSE, reduction = "tsne", split.by = "WA3")　　

DefaultAssay(WA3.combined.OL) <- "RNA"

#Marker top2
dataoligo.markers <- FindAllMarkers(dataoligo, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
dataoligo.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

#Heatmap
top5 <- dataoligo.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(dataoligo, features = top5$gene) + NoLegend()

#Top100
top100 <- dataoligo.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
top100pval <- subset(top100, rowSums(top100[5] < 0.05) > 0)

df <- top100pval[,7:6]
dfsample <- split(df$gene,df$cluster)
length(dfsample) 

BiocManager::install("AnnotationHub", force = TRUE)
library("AnnotationHub")

dfsample$`0` = bitr(dfsample$`0`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`1` = bitr(dfsample$`1`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`2` = bitr(dfsample$`2`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`3` = bitr(dfsample$`3`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`4` = bitr(dfsample$`4`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`5` = bitr(dfsample$`5`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

genelist <- list("0" = dfsample$`0`$ENTREZID,
                 "1" = dfsample$`1`$ENTREZID,
                 "2" = dfsample$`2`$ENTREZID,
                 "3" = dfsample$`3`$ENTREZID,
                 "4" = dfsample$`4`$ENTREZID,
                 "5" = dfsample$`5`$ENTREZID)
GOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Mm.eg.db")
dotplot(GOclusterplot)

#Reactome pathway
BiocManager::install("ReactomePA")
library(ReactomePA)
reactome <- compareCluster(geneClusters = genelist, 
                           fun = "enrichPathway",
                           organism ="mouse")
dotplot(reactome)  

######
##DEGs of each cell type
WA3.combined.subset <- readRDS("~/Downloads/R/230908 3w RPCA_before_annotation.rds")
WA3.combined.subset$cluster_sample <- paste(WA3.combined.subset$seurat_clusters, WA3.combined.subset$orig.ident, sep = "_")
Idents(WA3.combined.subset) <- "cluster_sample"
DEGcelltype_WA3 <- FindMarkers(WA3.combined.subset, ident.1 = c("XX_A3"), ident.2 = c("XX_W3"), logfc.threshold = 0.1)

