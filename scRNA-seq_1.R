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
W6.data <- Read10X_h5("W6-2nd-filtered_feature_bc_matrix.h5")
A6.data <- Read10X_h5("A6-2nd-filtered_feature_bc_matrix.h5")

#Creating the seurat object
W6 <- CreateSeuratObject(counts = W6.data,project = "W6", min.cell = 3)
A6 <- CreateSeuratObject(counts = A6.data,project = "A6", min.cell = 3)

#Store in a group named WA6
W6$WA6 <- "W6"
A6$WA6 <- "A6"

#QC and selecting cells for further analysis
W6[["percent.mt"]] <- PercentageFeatureSet(W6, pattern = "mt-")
W6 <- subset(W6, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 1)
A6[["percent.mt"]] <- PercentageFeatureSet(A6, pattern = "mt-")
A6 <- subset(A6, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 1)

#Make a list
WA6.list <- NULL
WA6.list <- list(W6, A6)
names(WA6.list) <- c("W6", "A6")

#Normalize and identify variable features for each dataset independently
WA6.list <- lapply(X = WA6.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000) 
})

#Select features that are repeatedly variable across datasets for integration run PCA on each
features <- SelectIntegrationFeatures(object.list = WA6.list)
WA6.list <- lapply(X = WA6.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

#Perform integration
WA6.anchors <- FindIntegrationAnchors(object.list = WA6.list, anchor.features = features, reduction = "rpca", k.anchor =2 )
WA6.combined <- IntegrateData(anchorset = WA6.anchors)
DefaultAssay(WA6.combined) <- "integrated"

#Visualization and clustering (tSNE)
WA6.combined <- ScaleData(WA6.combined, verbose = FALSE)
WA6.combined <- RunPCA(WA6.combined , npcs = 25, verbose = FALSE)
WA6.combined <- RunTSNE(WA6.combined , reduction = "pca", dims = 1:25)
WA6.combined  <- FindNeighbors(WA6.combined, graph.name = "WA6test", reduction = "pca", dims = 1:25)
WA6.combined  <- FindClusters(WA6.combined, graph.name = "WA6test", resolution = 0.8)

#Visualization and clustering (UMAP)
WA6.combined <- ScaleData(WA6.combined, verbose = FALSE)
WA6.combined <- RunPCA(WA6.combined , npcs = 25, verbose = FALSE)
WA6.combined <- RunUMAP(WA6.combined , reduction = "pca", dims = 1:25)
WA6.combined  <- FindNeighbors(WA6.combined, graph.name = "WA6test", reduction = "pca", dims = 1:25)
WA6.combined  <- FindClusters(WA6.combined, graph.name = "WA6test", resolution = 0.8)

DefaultAssay(WA6.combined) <- "RNA"  

#DotPlot 
features <- c("Pdgfra", "C1ql1", "Tnr", "Ust", "Cnksr3", "Ptgds", "Mobp", "Mog", "Mag", "Mbp", "Plp1", "Celf4", "Rbfox3", "Slc17a6", "Slc6a1", "Slc6a5", "Gad2","Slc5a7", "Chat","Slc1a2", "Aqp4", "Ly86", "Runx1", "Vtn", "Abcc9", "Flt1")
DotPlot(WA6.combined, features = features) 

#Exclude clusters containing doublets or background
WA6.combined.subset <- subset(WA6.combined,  idents = c("5"), invert=TRUE)
saveRDS(WA6.combined.subset,"/230907 6w RPCA_before_annotation.rds") 

#heatmap
WA6.combined.marker <- FindAllMarkers(WA6.combined.subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
WA6.combined.marker %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top5 <- WA6.combined.marker %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(WA6.combined.subset, features = top5$gene) + NoLegend()

##DotPlot
features <- c("Pdgfra", "C1ql1", "Tnr", "Ust", "Cnksr3", "Ptgds", "Mobp", "Mog", "Mag", "Mbp", "Plp1", "Celf4", "Rbfox3", "Slc17a6", "Slc6a1", "Slc6a5", "Gad2", "Chat","Slc1a2", "Aqp4", "Ly86", "Runx1", "Vtn", "Abcc9", "Flt1")
DotPlot(WA6.combined.subset, features = features) 

#Label clusters
new.cluster.ids <- c("Oligodendrocytes","Oligodendrocytes", "Astrocytes","Oligodendrocytes", 
                     "Oligodendrocytes","Inhibitory neurons", "OPC","Oligodendrocytes",
                     "Microglia","Inhibitory neurons", "Excitatory neurons","Excitatory neurons",
                     "Excitatory neurons","Astrocytes","Inhibitory neurons","Excitatory neurons",
                     "Excitatory neurons", "OL progenitors","Pericytes","Endothelial cells",
                     "Fibroblasts", "Inhibitory neurons")
names(new.cluster.ids) <- levels(WA6.combined.subset)
WA6.combined.subset <- RenameIdents(WA6.combined.subset, new.cluster.ids)
DimPlot(WA6.combined.subset , reduction = "tsne", label = TRUE, split.by="WA6")
DimPlot(WA6.combined.subset, reduction = "tsne", label = FALSE, 
        cols = c("Oligodendrocytes"="#F8766D","Inhibitory neurons" ="#CD9600", "Astrocytes"="#7CAE00", "Excitatory neurons"="#0CB702", "OPC"="#ED68ED", "OL progenitors"="#00A9FF", "Microglia"="#00BFC4", "Pericytes"="#8494FF", "Endothelial cells"="#C77CFF", "Fibroblasts"="#FF6C91" ),
        order = c("Fibroblasts","Endothelial cells","Pericytes","Microglia","OL progenitors","OPC","Astrocytes","Excitatory neurons", "Inhibitory neurons","Oligodendrocyte"))  ###採用Fig.3B
DimPlot(WA6.combined.subset, reduction = "umap", label = FALSE, 
        cols = c("Oligodendrocytes"="#F8766D","Inhibitory neurons" ="#CD9600", "Astrocytes"="#7CAE00", "Excitatory neurons"="#0CB702", "OPC"="#ED68ED", "OL progenitors"="#00A9FF", "Microglia"="#00BFC4", "Pericytes"="#8494FF", "Endothelial cells"="#C77CFF", "Fibroblasts"="#FF6C91" ),
        order = c("Fibroblasts","Endothelial cells","Pericytes","Microglia","OL progenitors","OPC","Astrocytes","Excitatory neurons", "Inhibitory neurons","Oligodendrocyte"))　###採用Fig.3B
table(Idents(WA6.combined.subset), WA6.combined.subset$WA6) 

#FeaturePlot
oligo <- subset(x = WA6.combined.subset, idents=c("Oligodendrocytes"))
features <- c("Asic2")  
FeaturePlot(object = oligo, features = features,split.by="WA6", reduction = "tsne",pt.size = 0.2,   min.cutoff = 1, max.cutoff = 3)
features <- c("Fam155a") 
FeaturePlot(object = oligo, features = features,split.by="WA6",reduction = "tsne", pt.size = 0.2,   min.cutoff = 1.5, max.cutoff =3)


WA6.combined.subset <- readRDS("230907 6w RPCA_before_annotation.rds")
#DEGs of OLs clusters
WA6.combined.subset$cluster_sample <- paste(WA6.combined.subset$seurat_clusters, WA6.combined.subset$orig.ident, sep = "_")
Idents(WA6.combined.subset) <- "cluster_sample"
DEGoligo_WA6 <- FindMarkers(WA6.combined.subset, ident.1 = c("0_A6","1_A6","3_A6","4_A6","8_A6"), ident.2 = c("0_W6","1_W6","3_W6","4_W6", "8_W6"), logfc.threshold = 0.1)

#GO enrichment analysis of 100 genes elevated in AR-97Q mice
DEGoligo_WA6_0.05 <- DEGoligo_WA6[DEGoligo_WA6$p_val_adj < 0.05
                                  & DEGoligo_WA6$avg_log2FC>
                                    0.404,]

DEG.gene_SYMBOLs <- rownames(DEGoligo_WA6_0.05)
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
DEGoligo_WA6_0.05_0.2 <- DEGoligo_WA6[DEGoligo_WA6$p_val_adj < 0.05
                                      & DEGoligo_WA6$avg_log2FC< -0.2237
                                      ,]

DEG.gene_SYMBOLs <- rownames(DEGoligo_WA6_0.05_0.2)
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
WA6.combined.subset <- readRDS("/230907 6w RPCA_before_annotation.rds")
DEGoligo_WA6 <- FindMarkers(WA6.combined.subset, ident.1 = c("0_A6","1_A6","3_A6","4_A6","8_A6"), ident.2 = c("0_W6","1_W6","3_W6","4_W6", "8_W6"), logfc.threshold = 0.1)

#Annotate the Ensembl gene IDs to gene symbols:
ens <- rownames(DEGoligo_WA6)
symbols <- mapIds(org.Mm.eg.db, keys = ens,
                  column = c('SYMBOL'), keytype = 'SYMBOL')
symbols <- symbols[!is.na(symbols)]
symbols <- symbols[match(rownames(DEGoligo_WA6), names(symbols))]
rownames(DEGoligo_WA6) <- symbols
keep <- !is.na(rownames(DEGoligo_WA6))
DEGoligo_WA6  <- DEGoligo_WA6 [keep,]

EnhancedVolcano(DEGoligo_WA6, 
                lab = rownames(DEGoligo_WA6),
                x ="avg_log2FC", 
                y ="p_val_adj",
                pCutoff = 10e-2, 
                FCcutoff = 0.2,
                pointSize = 1.0,selectLab = c( "Rbfox1","Meg3",'Asic2', "Kcnip4","Fam155a", "Usp31","Cdk19","Qk", "Car2", "4930419G24Rik"),
                labSize = 5.0,
                drawConnectors=TRUE, widthConnectors=0.5)   

#Reanalyse the OL cluster
WA6.combined.OL <- subset(WA6.combined.subset,  idents = c("0", "1", "3", "4", "8"))
all.genes <- rownames(WA6.combined.OL)
WA6.combined.OL <- ScaleData(WA6.combined.OL, features = all.genes)
WA6.combined.OL <- FindVariableFeatures(object = WA6.combined.OL)
WA6.combined.OL <- RunPCA(WA6.combined.OL, features = VariableFeatures(object = WA6.combined.OL))
WA6.combined.OL <- RunTSNE(WA6.combined.OL , reduction = "pca", dims = 1:25)
WA6.combined.OL  <- FindNeighbors(WA6.combined.OL, graph.name = "WA6test", reduction = "pca", dims = 1:25)
WA6.combined.OL  <- FindClusters(WA6.combined.OL, graph.name = "WA6test", resolution = 0.8)
p1 <- DimPlot(WA6.combined.OL , reduction = "tsne", group.by = "WA6")
plot_grid(p1)　

######
##DEGs of each cell type
WA6.combined.subset <- readRDS("/230907 6w RPCA_before_annotation.rds")
WA6.combined.subset$cluster_sample <- paste(WA6.combined.subset$seurat_clusters, WA6.combined.subset$orig.ident, sep = "_")
Idents(WA6.combined.subset) <- "cluster_sample"
DEGcelltype_WA6 <- FindMarkers(WA6.combined.subset, ident.1 = c("XX_A6"), ident.2 = c("XX_W6"), logfc.threshold = 0.1)

#Analysis of OPC clusters
WA6.combined.subset <- readRDS("230907 6w RPCA_before_annotation.rds")
WA6.combined.subset$cluster_sample <- paste(WA6.combined.subset$seurat_clusters, WA6.combined.subset$orig.ident, sep = "_")
Idents(WA6.combined.subset) <- "cluster_sample"
DEGopc_WA6 <- FindMarkers(WA6.combined.subset, ident.1 = c("7_A6"), ident.2 = c("7_W6"), logfc.threshold = 0.1)

#GO enrichment analysis of 19 genes upregulated in AR-97Q mice
DEGopc_WA6_0.05 <- DEGopc_WA6[DEGopc_WA6$p_val_adj < 0.05
                              & DEGopc_WA6$avg_log2FC>
                                0.3,]

DEG.gene_SYMBOLs <- rownames(DEGopc_WA6_0.05)
Mm <- org.Mm.eg.db
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
DEGopc_WA6 <- FindMarkers(WA6.combined.subset, ident.1 = c("7_A6"), ident.2 = c("7_W6"), logfc.threshold = 0.1)

#Annotate the Ensembl gene IDs to gene symbols:
ens <- rownames(DEGopc_WA6)
symbols <- mapIds(org.Mm.eg.db, keys = ens,
                  column = c('SYMBOL'), keytype = 'SYMBOL')
symbols <- symbols[!is.na(symbols)]
symbols <- symbols[match(rownames(DEGopc_WA6), names(symbols))]
rownames(DEGopc_WA6) <- symbols
keep <- !is.na(rownames(DEGopc_WA6))
DEGopc_WA6  <- DEGopc_WA6 [keep,]
EnhancedVolcano(DEGopc_WA6 , 
                lab = rownames(DEGopc_WA6),
                x ="avg_log2FC", 
                y ="p_val_adj",
                pCutoff = 10e-2, 
                FCcutoff = 0.2,
                pointSize = 1.0,selectLab = c("Kcnip4","Meg3",'Cntnap5a','Cntn5', "Kcnh7","Ppfibp1","Sema3d", "Nxph1", "Epn2","Fchsd2"),
                labSize = 5.0,
                drawConnectors=TRUE, widthConnectors=0.5)  

##Label clusters for Cellchat
new.cluster.ids <- c("OL","OL", "AS","OL",
                     "OL","IN", "OPC","OL",
                     "MI","IN", "EX","EX",
                     "EX","AS", "IN","EX",
                     "EX","OLpro", "PE",  "EN",
                     "FI", "IN")

names(new.cluster.ids) <- levels(WA6.combined.subset)
WA6.combined.subset <- RenameIdents(WA6.combined.subset, new.cluster.ids)
WA6.combined.subset$annotated <- Idents(WA6.combined.subset)
saveRDS(WA6.combined.subset,"~/Downloads/R/WA6.combined_anotated_231120.rds")

