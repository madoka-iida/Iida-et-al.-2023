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
W13.data <- Read10X_h5("W13 filtered_feature_bc_matrix-2.h5")
A13.data <- Read10X_h5("A13 filtered_feature_bc_matrix-2.h5")

#Creating the seurat object
W13 <- CreateSeuratObject(counts = W13.data,project = "W13", min.cell = 3)
A13 <- CreateSeuratObject(counts = A13.data,project = "A13", min.cell = 3)

#Store in a group named WA13
W13$WA13 <- "W13"
A13$WA13 <- "A13"

#QC and selecting cells for further analysis
W13[["percent.mt"]] <- PercentageFeatureSet(W13, pattern = "mt-")
W13 <- subset(W13, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 1)
A13[["percent.mt"]] <- PercentageFeatureSet(A13, pattern = "mt-")
A13<- subset(A13, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 1)

#Make a list
WA13.list <- NULL
WA13.list <- list(W13, A13)
names(WA13.list) <- c("W13", "A13")

#Normalize and identify variable features for each dataset independently
WA13.list <- lapply(X = WA13.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)  ##3000
}) 

#Select features that are repeatedly variable across datasets for integration run PCA on each
features <- SelectIntegrationFeatures(object.list = WA13.list)
WA13.list <- lapply(X = WA13.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

#Perform integration
WA13.anchors <- FindIntegrationAnchors(object.list = WA13.list, anchor.features = features, reduction = "rpca", k.anchor =2 )
WA13.combined <- IntegrateData(anchorset = WA13.anchors)
DefaultAssay(WA13.combined) <- "integrated"

#Visualization and clustering (tSNE)
WA13.combined <- ScaleData(WA13.combined, verbose = FALSE)
WA13.combined <- RunPCA(WA13.combined , npcs = 25, verbose = FALSE)
WA13.combined <- RunTSNE(WA13.combined , reduction = "pca", dims = 1:25)
WA13.combined  <- FindNeighbors(WA13.combined, graph.name = "WA13test", reduction = "pca", dims = 1:25)
WA13.combined  <- FindClusters(WA13.combined, graph.name = "WA13test", resolution = 0.8)

#Visualization and clustering (UMAP)
WA13.combined <- ScaleData(WA13.combined, verbose = FALSE)
WA13.combined <- RunPCA(WA13.combined , npcs = 25, verbose = FALSE)
WA13.combined <- RunUMAP(WA13.combined , reduction = "pca", dims = 1:25)
WA13.combined  <- FindNeighbors(WA13.combined, graph.name = "WA13test", reduction = "umap", dims = 1:25)
WA13.combined  <- FindClusters(WA13.combined, graph.name = "WA13test", resolution = 0.8)

DefaultAssay(WA13.combined) <- "RNA" 

#DotPlot
features <- c("Pdgfra", "C1ql1", "Tnr", "Ust", "Cnksr3", "Ptgds", "Mobp", "Mog", "Mag", "Mbp", "Plp1", "Celf4", "Rbfox3", "Slc17a6", "Slc6a1", "Slc6a5", "Gad2","Slc5a7", "Chat","Slc1a2", "Aqp4", "Ly86", "Runx1", "Vtn", "Abcc9", "Flt1")
DotPlot(WA13.combined, features = features) 

#Exclude clusters containing doublets or background
WA13.combined.subset <- subset(WA13.combined,  idents = c("6", "22"), invert=TRUE)
saveRDS(WA13.combined.subset,"/230907 13w RPCA_before_annotation.rds") 

#heatmap
WA13.combined.marker <- FindAllMarkers(WA13.combined.subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
WA13.combined.marker %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top5 <- WA13.combined.marker %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(WA13.combined.subset, features = top5$gene) + NoLegend()

#DotPlot 
features <- c("Pdgfra", "C1ql1", "Tnr", "Ust", "Cnksr3", "Ptgds", "Mobp", "Mog", "Mag", "Mbp", "Plp1", "Celf4", "Rbfox3", "Slc17a6", "Slc6a1", "Slc6a5", "Gad2", "Chat","Slc1a2", "Aqp4", "Ly86", "Runx1", "Vtn", "Abcc9", "Flt1")
DotPlot(WA13.combined.subset, features = features) 


#Label clusters 
new.cluster.ids <- c("Oligodendrocytes","Oligodendrocytes","Astrocytes", "Oligodendrocytes", 
                     "Oligodendrocytes", "Oligodendrocytes", "OPC","Inhibitory neurons",
                     "Microglia","Excitatory neurons","Microglia","Excitatory neurons", 
                     "Oligodendrocytes","Excitatory neurons","Astrocytes","OL progenitors", 
                     "Excitatory neurons", "Excitatory neurons","Pericytes", "Inhibitory neurons", 
                     "Inhibitory neurons", "Fibroblasts")
names(new.cluster.ids) <- levels(WA13.combined.subset)
WA13.combined.subset <- RenameIdents(WA13.combined.subset, new.cluster.ids)
DimPlot(WA13.combined.subset , reduction = "tsne", label = TRUE, split.by="WA13")
DimPlot(WA13.combined.subset, reduction = "tsne", label = FALSE,
        cols = c("Oligodendrocytes"="#F8766D","Inhibitory neurons" ="#CD9600", "Astrocytes"="#7CAE00", "Excitatory neurons"="#0CB702", "OPC"="#ED68ED", "OL progenitors"="#00A9FF", "Microglia"="#00BFC4", "Pericytes"="#8494FF", "Fibroblasts"="#FF6C91"),
        order = c("Fibroblasts","Pericytes","Microglia","OL progenitors","OPC","Astrocytes","Excitatory neurons", "Inhibitory neurons","Oligodendrocyte"))  
DimPlot(WA13.combined.subset, reduction = "umap", label = FALSE,
        cols = c("Oligodendrocytes"="#F8766D","Inhibitory neurons" ="#CD9600", "Astrocytes"="#7CAE00", "Excitatory neurons"="#0CB702", "OPC"="#ED68ED", "OL progenitors"="#00A9FF", "Microglia"="#00BFC4", "Pericytes"="#8494FF", "Fibroblasts"="#FF6C91"),
        order = c("Fibroblasts","Pericytes","Microglia","OL progenitors","OPC","Astrocytes","Excitatory neurons", "Inhibitory neurons","Oligodendrocyte")) 
table(Idents(WA13.combined.subset), WA13.combined.subset$WA13) 

#FeaturePlot
oligo <- subset(x = WA13.combined.subset, idents=c("Oligodendrocytes"))
features <- c("Asic2") 
FeaturePlot(object = oligo, features = features,reduction = "tsne",split.by="WA13", pt.size = 0.2,   min.cutoff = 1, max.cutoff = 3)
features <- c("Fam155a") 
FeaturePlot(object = oligo, features = features,reduction = "tsne",split.by="WA13", pt.size = 0.2,   min.cutoff = 1.5, max.cutoff =3.5)


WA13.combined.subset <- readRDS("/230907 13w RPCA_before_annotation.rds")
#DEGs of OLs clusters
WA13.combined.subset$cluster_sample <- paste(WA13.combined.subset$seurat_clusters, WA13.combined.subset$orig.ident, sep = "_")
Idents(WA13.combined.subset) <- "cluster_sample"
DEGoligo_WA13 <- FindMarkers(WA13.combined.subset, ident.1 = c("0_A13","1_A13","3_A13","4_A13","5_A13","13_A13"), ident.2 = c("0_W13","1_W13","3_W13","4_W13","5_W13","13_W13"), logfc.threshold = 0.1)

#GO enrichment analysis of 100 genes elevated in AR-97Q mice
DEGoligo_WA13_0.05 <- DEGoligo_WA13[DEGoligo_WA13$p_val_adj < 0.05
                                    & DEGoligo_WA13$avg_log2FC>
                                      0.415,]
DEG.gene_SYMBOLs <- rownames(DEGoligo_WA13_0.05)
Mm <- org.Mm.eg.db
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


#GO enrichment analysis of 100 genes downregulated in AR-97Q mice
DEGoligo_WA13_0.05 <- DEGoligo_WA13[DEGoligo_WA13$p_val_adj < 0.05
                                    & DEGoligo_WA13$avg_log2FC< -0.436
                                    ,]

DEG.gene_SYMBOLs <- rownames(DEGoligo_WA13_0.05)
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
WA13.combined.subset <- readRDS("/230907 13w RPCA_before_annotation.rds")
DEGoligo_WA13 <- FindMarkers(WA13.combined.subset, ident.1 = c("0_A13","1_A13","3_A13","4_A13","5_A13","13_A13"), ident.2 = c("0_W13","1_W13","3_W13","4_W13","5_W13","13_W13"), logfc.threshold = 0.1)

#Annotate the Ensembl gene IDs to gene symbols:
ens <- rownames(DEGoligo_WA13)
symbols <- mapIds(org.Mm.eg.db, keys = ens,
                  column = c('SYMBOL'), keytype = 'SYMBOL')
symbols <- symbols[!is.na(symbols)]
symbols <- symbols[match(rownames(DEGoligo_WA13), names(symbols))]
rownames(DEGoligo_WA13) <- symbols
keep <- !is.na(rownames(DEGoligo_WA13))
DEGoligo_WA13  <- DEGoligo_WA13 [keep,]

EnhancedVolcano(DEGoligo_WA13, 
                lab = rownames(DEGoligo_WA13),
                x ="avg_log2FC", 
                y ="p_val_adj",
                pCutoff = 10e-2, 
                FCcutoff = 0.2,
                pointSize = 1.0,selectLab = c( "9330111N05Rik","Garnl3",'A330015K06Rik', "Kcnma1","Rhoj", "Dpp10","Sgcz","Nrg3", "Cntnap2", "P4ha1" ),
                labSize = 5.0,
                drawConnectors=TRUE, widthConnectors=0.5)   

#Reanalyze the OL cluster
WA13.combined.OL <- subset(WA13.combined.subset,  idents = c("0", "1", "3","4", "5", "13"))
all.genes <- rownames(WA13.combined.OL)
WA13.combined.OL <- ScaleData(WA13.combined.OL, features = all.genes)
WA13.combined.OL <- FindVariableFeatures(object = WA13.combined.OL)
WA13.combined.OL <- RunPCA(WA13.combined.OL, features = VariableFeatures(object = WA13.combined.OL))
WA13.combined.OL <- RunTSNE(WA13.combined.OL , reduction = "pca", dims = 1:25)
WA13.combined.OL  <- FindNeighbors(WA13.combined.OL, graph.name = "WA13test", reduction = "pca", dims = 1:25)
WA13.combined.OL  <- FindClusters(WA13.combined.OL, graph.name = "WA13test", resolution = 0.8)
p1 <- DimPlot(WA13.combined.OL , reduction = "tsne", group.by = "WA13")
plot_grid(p1)　

######
##DEGs of each cell type
WA13.combined.subset <- readRDS("/230907 13w RPCA_before_annotation.rds")
WA13.combined.subset$cluster_sample <- paste(WA13.combined.subset$seurat_clusters, WA13.combined.subset$orig.ident, sep = "_")
Idents(WA13.combined.subset) <- "cluster_sample"
DEGcelltype_WA13 <- FindMarkers(WA13.combined.subset, ident.1 = c("XX_A13"), ident.2 = c("XX_W13"), logfc.threshold = 0.1)

#Analysis of OPC clusters
WA13.combined.subset <- readRDS("~/Downloads/R/230907 13w RPCA_before_annotation.rds")
WA13.combined.subset$cluster_sample <- paste(WA13.combined.subset$seurat_clusters, WA13.combined.subset$orig.ident, sep = "_")
Idents(WA13.combined.subset) <- "cluster_sample"
DEGopc_WA13 <- FindMarkers(WA13.combined.subset, ident.1 = c("7_A13"), ident.2 = c("7_W13"), logfc.threshold = 0.1)

#GO enrichment analysis of 100 genes suppressed in AR-97Q mice
DEGopc_WA13_0.05 <- DEGopc_WA13[DEGopc_WA13$p_val_adj < 0.05
                                & DEGopc_WA13$avg_log2FC< -0.372
                                ,]
DEG.gene_SYMBOLs <- rownames(DEGopc_WA13_0.05)
Mm<- org.Mm.eg.db
DEG.gene_IDs <- AnnotationDbi::select(Mm,
                                      keys=DEG.gene_SYMBOLs,
                                      columns = c(
                                        "SYMBOL",
                                        "ENTREZID"
                                      ),
                                      keytype= "SYMBOL"  )$ENTREZID
#BP
ego <- enrichGO(gene = DEG.gene_IDs,
                ont ="BP",
                OrgDb='org.Mm.eg.db')
dotplot(ego, showCategory=5)　

#Volcano plot of OPC clusters
DEGopc_WA13 <- FindMarkers(WA13.combined.subset, ident.1 = c("7_A13"), ident.2 = c("7_W13"), logfc.threshold = 0.1)

# Annotate the Ensembl gene IDs to gene symbols:
ens <- rownames(DEGopc_WA13)

symbols <- mapIds(org.Mm.eg.db, keys = ens,
                  column = c('SYMBOL'), keytype = 'SYMBOL')

symbols <- symbols[!is.na(symbols)]
symbols <- symbols[match(rownames(DEGopc_WA13), names(symbols))]
rownames(DEGopc_WA13) <- symbols
keep <- !is.na(rownames(DEGopc_WA13))
DEGopc_WA13  <- DEGopc_WA13 [keep,]
EnhancedVolcano(DEGopc_WA13, 
                lab = rownames(DEGopc_WA13),
                x ="avg_log2FC", 
                y ="p_val_adj",
                pCutoff = 10e-2,  ###自由に変更する
                FCcutoff = 0.2,
                pointSize = 1.0,selectLab = c( "Kcnip1",'Mthfd1l', "Neat1","Spata13", "9330111N05Rik", "Lingo2","Cntnap2","Gm47283", "Cntn5" ),
                labSize = 5.0,　　#6位にkcnip4
                drawConnectors=TRUE, widthConnectors=0.5) 

##Label clusters　for CellChat
new.cluster.ids <- c("OL","OL","AS", "OL", 
                     "OL","OL", "OPC","IN",
                     "MI","EX","MI", "EX",
                     "OL", "EX","AS","OLpro", 
                     "EX", "EX", "PE", "IN",
                     "IN", "FI")
