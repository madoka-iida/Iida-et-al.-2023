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
W9.data <- Read10X_h5("W9 filtered_feature_bc_matrix-2.h5")
A9.data <- Read10X_h5("A9 filtered_feature_bc_matrix-2.h5")

#Creating the seurat object
W9 <- CreateSeuratObject(counts = W9.data,project = "W9", min.cell = 3)
A9 <- CreateSeuratObject(counts = A9.data,project = "A9", min.cell = 3)

#Store in a group named WA9
W9$WA9 <- "W9"
A9$WA9 <- "A9"

#QC and selecting cells for further analysis
W9[["percent.mt"]] <- PercentageFeatureSet(W9, pattern = "mt-")
W9 <- subset(W9, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 1)
A9[["percent.mt"]] <- PercentageFeatureSet(A9, pattern = "mt-")
A9 <- subset(A9, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 1)

#Make a list
WA9.list <- NULL
WA9.list <- list(W9, A9)
names(WA9.list) <- c("W9", "A9")

#Normalize and identify variable features for each dataset independently
WA9.list <- lapply(X = WA9.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)  #### 2000ではなく3000
})

#Select features that are repeatedly variable across datasets for integration run PCA on each
features <- SelectIntegrationFeatures(object.list = WA9.list)
WA9.list <- lapply(X = WA9.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

#Perform integration
WA9.anchors <- FindIntegrationAnchors(object.list = WA9.list, anchor.features = features, reduction = "rpca", k.anchor =2 )
WA9.combined <- IntegrateData(anchorset = WA9.anchors)
DefaultAssay(WA9.combined) <- "integrated"

#Visualization and clustering (tSNE)
WA9.combined <- ScaleData(WA9.combined, verbose = FALSE)
WA9.combined <- RunPCA(WA9.combined , npcs = 25, verbose = FALSE)
WA9.combined <- RunTSNE(WA9.combined , reduction = "pca", dims = 1:25)
WA9.combined  <- FindNeighbors(WA9.combined, graph.name = "WA9test", reduction = "pca", dims = 1:25)
WA9.combined  <- FindClusters(WA9.combined, graph.name = "WA9test", resolution = 0.8)

#Visualization and clustering (UMAP)
WA9.combined <- ScaleData(WA9.combined, verbose = FALSE)
WA9.combined <- RunPCA(WA9.combined , npcs = 25, verbose = FALSE)
WA9.combined <- RunUMAP(WA9.combined , reduction = "pca", dims = 1:25)
WA9.combined  <- FindNeighbors(WA9.combined, graph.name = "WA9test", reduction = "pca", dims = 1:25)
WA9.combined  <- FindClusters(WA9.combined, graph.name = "WA9test", resolution = 0.8)

DefaultAssay(WA9.combined) <- "RNA"  

#DotPlot
features <- c("Pdgfra", "C1ql1", "Tnr", "Ust", "Cnksr3", "Ptgds", "Mobp", "Mog", "Mag", "Mbp", "Plp1", "Celf4", "Rbfox3", "Slc17a6", "Slc6a1", "Slc6a5", "Gad2", "Chat","Slc1a2", "Aqp4", "Ly86", "Runx1", "Vtn", "Abcc9", "Flt1","Dnah12")
DotPlot(WA9.combined, features = features) 

#Exclude clusters containing doublets or background
WA9.combined.subset <- subset(WA9.combined,  idents = c("6"), invert=TRUE)
saveRDS(WA9.combined.subset,"~/Downloads/R/230907 9w RPCA_before_annotation.rds") 

#Heatmap
WA9.combined.marker <- FindAllMarkers(WA9.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
WA9.combinedmarker %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top5 <- WA9.combined.marker %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(WA9.combined, features = top5$gene) + NoLegend()

#DotPlot 
features <- c("Pdgfra", "C1ql1", "Tnr", "Ust", "Cnksr3", "Ptgds", "Mobp", "Mog", "Mag", "Mbp", "Plp1", "Celf4", "Rbfox3", "Slc17a6", "Slc6a1", "Slc6a5", "Gad2", "Chat","Slc1a2", "Aqp4", "Ly86", "Runx1", "Vtn", "Abcc9", "Flt1")
DotPlot(WA9.combined.subset, features = features) 


#Label clusters 
new.cluster.ids <- c("Oligodendrocytes","Oligodendrocytes", "Inhibitory neurons","Astrocytes",
                     "Oligodendrocytes","Oligodendrocytes",  "OPC","Microglia", 
                     "Inhibitory neurons","Excitatory neurons", "Excitatory neurons", "Inhibitory neurons",
                     "Excitatory neurons","Oligodendrocytes","Excitatory neurons","Excitatory neurons",
                     "Excitatory neurons","OL progenitors","Excitatory neurons","Pericytes", 
                     "Inhibitory neurons","Endothelial cells")
names(new.cluster.ids) <- levels(WA9.combined.subset)
WA9.combined.subset <- RenameIdents(WA9.combined.subset, new.cluster.ids)
DimPlot(WA9.combined.subset , reduction = "tsne", label = TRUE, split.by="WA9")
DimPlot(WA9.combined.subset, reduction = "tsne", label = FALSE, 
        cols = c("Oligodendrocytes"="#F8766D","Inhibitory neurons" ="#CD9600", "Astrocytes"="#7CAE00", "Excitatory neurons"="#0CB702", "OPC"="#ED68ED", "OL progenitors"="#00A9FF", "Microglia"="#00BFC4", "Pericytes"="#8494FF", "Endothelial cells"="#C77CFF" ),
        order = c("Endothelial cells","Pericytes","Microglia","OL progenitors","OPC","Astrocytes","Excitatory neurons", "Inhibitory neurons","Oligodendrocyte"))   
DimPlot(WA9.combined.subset, reduction = "umap", label = FALSE, 
        cols = c("Oligodendrocytes"="#F8766D","Inhibitory neurons" ="#CD9600", "Astrocytes"="#7CAE00", "Excitatory neurons"="#0CB702", "OPC"="#ED68ED", "OL progenitors"="#00A9FF", "Microglia"="#00BFC4", "Pericytes"="#8494FF", "Endothelial cells"="#C77CFF" ),
        order = c("Endothelial cells","Pericytes","Microglia","OL progenitors","OPC","Astrocytes","Excitatory neurons", "Inhibitory neurons","Oligodendrocyte"))    
table(Idents(WA9.combined.subset), WA9.combined.subset$WA9)  

#FeaturePlot
oligo <- subset(x = WA9.combined.subset, idents=c("Oligodendrocytes"))
features <- c("Asic2")  
FeaturePlot(object = oligo, features = features,reduction = "tsne",split.by="WA9", pt.size = 0.2,   min.cutoff = 1, max.cutoff = 4)
features <- c("Fam155a")　
FeaturePlot(object = oligo, features = features,reduction = "tsne",split.by="WA9", pt.size = 0.2,   min.cutoff = 2, max.cutoff =4)


WA9.combined.subset <- readRDS("~/Downloads/R/230907 9w RPCA_before_annotation.rds")
#DEGs of OLs clusters
WA9.combined$cluster_sample <- paste(WA9.combined$seurat_clusters, WA9.combined$orig.ident, sep = "_")
Idents(WA9.combined) <- "cluster_sample"
DEGoligo_WA9 <- FindMarkers(WA9.combined, ident.1 = c("0_A9","1_A9","4_A9","5_A9","14_A9"), ident.2 = c("0_W9","1_W9","4_W9","5_W9","14_W9"),logfc.threshold = 0.1)

#GO enrichment analysis of 100 genes elevated in AR-97Q mice
DEGoligo_WA9_0.05 <- DEGoligo_WA9[DEGoligo_WA9$p_val_adj < 0.05
                                  & DEGoligo_WA9$avg_log2FC>
                                    0.338,]
Mm <- org.Mm.eg.db
DEG.gene_SYMBOLs <- rownames(DEGoligo_WA9_0.05)
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
DEGoligo_WA9_0.05_0.2 <- DEGoligo_WA9[DEGoligo_WA9$p_val_adj < 0.05
                                      & DEGoligo_WA9$avg_log2FC< -0.234
                                      ,]

DEG.gene_SYMBOLs <- rownames(DEGoligo_WA9_0.05_0.2)
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
WA9.combined.subset <- readRDS("~/Downloads/R/230907 9w RPCA_before_annotation.rds")
DEGoligo_WA9 <- FindMarkers(WA9.combined, ident.1 = c("0_A9","1_A9","4_A9","5_A9","14_A9"), ident.2 = c("0_W9","1_W9","4_W9","5_W9","14_W9"),logfc.threshold = 0.1)

#Annotate the Ensembl gene IDs to gene symbols:
ens <- rownames(DEGoligo_WA9)
symbols <- mapIds(org.Mm.eg.db, keys = ens,
                  column = c('SYMBOL'), keytype = 'SYMBOL')
symbols <- symbols[!is.na(symbols)]
symbols <- symbols[match(rownames(DEGoligo_WA9), names(symbols))]
rownames(DEGoligo_WA9) <- symbols
keep <- !is.na(rownames(DEGoligo_WA9))
DEGoligo_WA9  <- DEGoligo_WA9 [keep,]

EnhancedVolcano(DEGoligo_WA9 , 
                lab = rownames(DEGoligo_WA9),
                x ="avg_log2FC", 
                y ="p_val_adj",
                pCutoff = 10e-2,  
                FCcutoff = 0.2,
                pointSize = 1.0,selectLab = c( "Asic2",'Meg3',"Etl4", "Rbfox1","Fam155a", "Fth1","Nrbp2", "Cldn11","Opalin", "St8sia3os"),
                labSize = 5.0,
                drawConnectors=TRUE, widthConnectors=0.5)  

#Reanalyze the OL cluster
WA9.combined.OL <- subset(WA9.combined.subset,  idents = c("0", "1", "4", "5", "14"))
all.genes <- rownames(WA9.combined.OL)
WA9.combined.OL <- ScaleData(WA9.combined.OL, features = all.genes)
WA9.combined.OL <- FindVariableFeatures(object = WA9.combined.OL)
WA9.combined.OL <- RunPCA(WA9.combined.OL, features = VariableFeatures(object = WA9.combined.OL))
WA9.combined.OL <- RunTSNE(WA9.combined.OL , reduction = "pca", dims = 1:25)
WA9.combined.OL  <- FindNeighbors(WA9.combined.OL, graph.name = "WA9test", reduction = "pca", dims = 1:25)
WA9.combined.OL  <- FindClusters(WA9.combined.OL, graph.name = "WA9test", resolution = 0.8)
p1 <- DimPlot(WA9.combined.OL , reduction = "tsne", group.by = "WA9")
plot_grid(p1)　

######
##DEGs of each cell type
WA9.combined.subset <- readRDS("~/Downloads/R/230907 9w RPCA_before_annotation.rds")
WA9.combined.subset$cluster_sample <- paste(WA9.combined.subset$seurat_clusters, WA9.combined.subset$orig.ident, sep = "_")
Idents(WA9.combined.subset) <- "cluster_sample"
DEGcelltype_WA9 <- FindMarkers(WA9.combined.subset, ident.1 = c("XX_A9"), ident.2 = c("XX_W9"), logfc.threshold = 0.1)

##Label clusters for Cellchat
new.cluster.ids <- c("OL","OL", "IN","AS",
                     "OL","OL", "OPC","MI",
                     "IN", "EX", "EX", "IN",
                     "EX","OL","EX","EX",
                     "EX","OLpro","EX","PE",
                     "IN","EN")

names(new.cluster.ids) <- levels(WA9.combined.subset)
WA9.combined.subset <- RenameIdents(WA9.combined.subset, new.cluster.ids)
WA9.combined.subset$annotated <- Idents(WA9.combined.subset)
saveRDS(WA9.combined.subset,"~/Downloads/R/WA9.combined_anotated_231120.rds")


