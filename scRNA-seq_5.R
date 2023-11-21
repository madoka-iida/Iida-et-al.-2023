#The scRNA-seq dataset was processed and analyzed with Seurat version 3.1.460.
# Load libraries
library(Seurat) #version 3.1.460
library(cowplot)
library(dplyr)
library(ggplot2)

#Loading the dataset
W_3.data <- Read10X_h5("W3 filtered_feature_bc_matrix-2.h5")
A_3.data <- Read10X_h5("A3 filtered_feature_bc_matrix-2.h5") #Other data as well.

#Creating the seurat object
W_3 <- CreateSeuratObject(counts = W_3.data,project = "W_3", min.cell = 3)
A_3 <- CreateSeuratObject(counts = A_3.data,project = "A_3", min.cell = 3)　#Other data as well.

#Store in a group named All
W_3$All <- "W_3"
A_3$All <- "A_3"　#Other data as well.

##QC and selecting cells for further analysis
W_3[["percent.mt"]] <- PercentageFeatureSet(W_3, pattern = "mt-")
W_3 <- subset(W_3, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 1)
A_3[["percent.mt"]] <- PercentageFeatureSet(A_3, pattern = "mt-")
A_3 <- subset(A_3, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 1) #Other data as well.

#Visualization
VlnPlot(W3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.05)
VlnPlot(A3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.05) #Other data as well.

#Normalization and detection of variable genes
W_3 <- NormalizeData(W_3, normalization.method = "LogNormalize", scale.factor = 10000)
A_3 <- NormalizeData(A_3, normalization.method = "LogNormalize", scale.factor = 10000)
W_3 <- FindVariableFeatures(W_3, selection.method = "vst", nfeatures = 3000)
A_3 <- FindVariableFeatures(A_3, selection.method = "vst", nfeatures = 3000) #Other data as well.

#Make lists
WA1.list <- NULL
WA2.list <- NULL
WA3.list <- NULL
WA4.list <- NULL

WA1.list <- list(W_3, A_3)
WA2.list <- list(W_6, A_6)
WA3.list <- list(W_9, A_9)
WA4.list <- list(W_13, A_13)

names(WA1.list) <- c("W_3", "A_3")
names(WA2.list) <- c("W_6", "A_6")
names(WA3.list) <- c("W_9", "A_9")
names(WA4.list) <- c("W_13", "A_13")

##Integrate normalized dataset
W_3 <- NormalizeData(W_3)
A_3 <- NormalizeData(A_3)
WA1.integrated <- merge(W_3, y = A_3, add.cell.ids = c('W_3', 'A_3'), project = 'WA1', merge.data = TRUE) #Other data as well.

##All_2にセット
WA1.integrated$All_2 <- "WA1"
WA2.integrated$All_2 <- "WA2"
WA3.integrated$All_2 <- "WA3"
WA4.integrated$All_2 <- "WA4"

##make a list
All_2.list <- NULL
All_2.list <- list(WA1.integrated, WA2.integrated, WA3.integrated, WA4.integrated)
names(All_2.list) <- c("WA1", "WA2", "WA3", "WA4")

#Integrate dataset
All_2.integrated <- merge(WA1.integrated, y = c(WA2.integrated, WA3.integrated, WA4.integrated), add.cell.ids = c("WA1", "WA2", "WA3", "WA4"), project = 'All_2', merge.data = TRUE)

##PCA
all_2.genes <- rownames(All_2.integrated)
All_2.integrated <- ScaleData(All_2.integrated, features = all_2.genes)
All_2.integrated <- FindVariableFeatures(object = All_2.integrated)
All_2.integrated <- RunPCA(All_2.integrated, features = VariableFeatures(object = All_2.integrated))

#Run Non-linear dimensional reduction (tSNE)
All_2.integrated <- ScaleData(All_2.integrated, verbose = FALSE)
All_2.integrated <- RunPCA(All_2.integrated , npcs = 25, verbose = FALSE)
All_2.integrated <- RunTSNE(All_2.integrated , reduction = "pca", dims = 1:25)
All_2.integrated  <- FindNeighbors(All_2.integrated, graph.name = "All2test", reduction = "pca", dims = 1:25)
All_2.integrated  <- FindClusters(All_2.integrated, graph.name = "All2test", resolution = 0.8)

#Run Non-linear dimensional reduction (UMAP)
All_2.integrated <- ScaleData(All_2.integrated, verbose = FALSE)
All_2.integrated <- RunPCA(All_2.integrated , npcs = 25, verbose = FALSE)
All_2.integrated <- RunUMAP(All_2.integrated , reduction = "pca", dims = 1:25)
All_2.integrated  <- FindNeighbors(All_2.integrated, graph.name = "All2test", reduction = "pca", dims = 1:25)
All_2.integrated  <- FindClusters(All_2.integrated, graph.name = "All2test", resolution = 0.8)

#Cluster identification  
All_2.integrated.marker <- FindAllMarkers(All_2.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
All_2.integrated.marker %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top10 <- All_2.integrated.marker %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(All_2.integrated, features = top10$gene) + NoLegend()　　

DefaultAssay(All_2.integrated) <- "RNA" 

#DotPlot
features <- c("Pdgfra", "C1ql1", "Tnr", "Ust", "Cnksr3", "Ptgds", "Mobp", "Mog", "Mag", "Mbp", "Plp1", "Celf4", "Rbfox3", "Slc17a6", "Slc6a1", "Slc6a5", "Gad2","Slc5a7", "Chat","Slc1a2", "Aqp4", "Ly86", "Runx1", "Vtn", "Abcc9", "Flt1", "Dnah12")
DotPlot(All_2.integrated, features = features) 　

#Exclude clusters containing doublets or background
All_2.integrated.subset <- subset(All_2.integrated,  idents = c("8","21", "22", "23"), invert=TRUE)
p1 <- DimPlot(All_2.integrated.subset, reduction = "tsne",  pt.size = 0.1, group.by = "All")
plot_grid(p1)　　

#Label clusters
new.cluster.ids <- c("Oligodendrocytes","Oligodendrocytes","Inhibitory neurons","Oligodendrocytes", #0-3
                     "Oligodendrocytes","Astrocytes","Oligodendrocytes","OPC", #4-7
                     "Microglia", "Excitatory neurons", "Inhibitory neurons","Excitatory neurons", #9-12
                     "Excitatory neurons","Oligodendrocytes","Excitatory neurons","OL progenitors", #13-16
                     "Inhibitory neurons","Astrocytes","Excitatory neurons","Pericytes",  #17-20
                     "Endothelial cells","Inhibitory neurons", "Fibroblasts", "Excitatory neurons") #24-27
names(new.cluster.ids) <- levels(All_2.integrated.subset)
All_2.integrated.subset <- RenameIdents(All_2.integrated.subset, new.cluster.ids)

DimPlot(All_2.integrated.subset, reduction = "tsne", 　
        cols = c("Oligodendrocytes"="#F8766D","Inhibitory neurons" ="#CD9600", "Astrocytes"="#7CAE00", "Excitatory neurons"="#0CB702", "OPC"="#ED68ED", "OL progenitors"="#00A9FF", "Microglia"="#00BFC4", "Pericytes"="#8494FF", "Endothelial cells"="#C77CFF", "Fibroblasts"="#FF6C91" ),
        order = c("Fibroblasts","Endothelial cells","Pericytes","Microglia","OL progenitors","OPC","Astrocytes","Excitatory neurons", "Inhibitory neurons","Oligodendrocyte"))　　####採用Fig.1B
DimPlot(All_2.integrated.subset, reduction = "umap",  　
        cols = c("Oligodendrocytes"="#F8766D","Inhibitory neurons" ="#CD9600", "Astrocytes"="#7CAE00", "Excitatory neurons"="#0CB702", "OPC"="#ED68ED", "OL progenitors"="#00A9FF", "Microglia"="#00BFC4", "Pericytes"="#8494FF", "Endothelial cells"="#C77CFF", "Fibroblasts"="#FF6C91" ),
        order = c("Fibroblasts","Endothelial cells","Pericytes","Microglia","OL progenitors","OPC","Astrocytes","Excitatory neurons", "Inhibitory neurons","Oligodendrocyte"))　　####採用Fig.1B

table(Idents(All_2.integrated.subset), All_2.integrated.subset$All)  

##DotPlot
levels(All_2.integrated.subset) <- c("Endothelial cells","Pericytes","Microglia","OL progenitors","OPC","Astrocytes","Excitatory neurons", "Inhibitory neurons","Oligodendrocytes")
DotPlot(All_2.integrated.subset, features = features) 　

#OL differentiation markers
features <- c("Pdgfra","Cspg4", "Neu4", "Sox6", "Bmp4", "Gpr17")
levels(All_2.integrated.subset) <- c("OPC","OL progenitors","Oligodendrocytes","","","","","","")
DotPlot(All_2.integrated.subset, features = features, idents =c("OPC","OL progenitors","Oligodendrocytes") ) 　####採用Fig.1C(Suppl?)

