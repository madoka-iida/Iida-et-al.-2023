#The scRNA-seq dataset was processed and analyzed with Seurat version 3.1.460.
#Load libraries
require(devtools)
install_version("spatstat.utils", version = "3.0.2")
library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(rvcheck)
library(org.Mm.eg.db)
library(EnhancedVolcano)
library(gridExtra)
library(grid)

#Loading the dataset, QC, and selecting cells for further analysis
AA.data <- Read10X_h5("A3 filtered_feature_bc_matrix-2.h5")
AB.data <- Read10X_h5("A6-2nd-filtered_feature_bc_matrix.h5")
AC.data <- Read10X_h5("A9 filtered_feature_bc_matrix-2.h5")
AD.data <- Read10X_h5("A13 filtered_feature_bc_matrix-2.h5")

#QC and selecting cells for further analysis

#Normalization and detection of variable genes
AA <- NormalizeData(AA, normalization.method = "LogNormalize", scale.factor = 10000)　#Other data as well.
AA <- FindVariableFeatures(AA, selection.method = "vst", nfeatures = 3000) #Other data as well.

#Make lists
AR.list <- NULL
AR.list <- list(AA, AB, AC, AD)
names(AR.list) <- c("AA", "AB", "AC", "AD")

##Integrate normalized dataset
AA <- NormalizeData(AA)　#Other data as well.
AR.integrated <- merge(AA, y = c(AB, AC, AD), add.cell.ids = c('AA', 'AB', 'AC', 'AD'), project = 'AR', merge.data = TRUE)
GetAssayData(AR.integrated)[1:10, 1:15]

#save
saveRDS(AR.integrated,"~/Downloads/R/230908 A36913 merge DEG_after_merge.rds")   #For subclustering
AR.integrated <- readRDS("~/Downloads/R/230908 A36913 merge DEG_after_merge.rds")

##PCA
all.genes <- rownames(AR.integrated)
AR.integrated <- ScaleData(AR.integrated, features = all.genes)
AR.integrated <- FindVariableFeatures(object = AR.integrated)
AR.integrated <- RunPCA(AR.integrated, features = VariableFeatures(object = AR.integrated))

#Run Non-linear dimensional reduction (UMAP)
AR.integrated <- ScaleData(AR.integrated, verbose = FALSE)
AR.integrated <- RunPCA(AR.integrated , npcs = 25, verbose = FALSE)
AR.integrated <- RunUMAP(AR.integrated , reduction = "pca", dims = 1:25)
AR.integrated  <- FindNeighbors(AR.integrated, graph.name = "ARtest", reduction = "pca", dims = 1:25)
AR.integrated  <- FindClusters(AR.integrated, graph.name = "ARtest", resolution = 1.2)

##Cluster identification
AR.integrated.marker <- FindAllMarkers(AR.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
AR.integrated.marker %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top10 <- AR.integrated.marker %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(AR.integrated, features = top10$gene) + NoLegend()

DefaultAssay(AR.integrated) <- "RNA" 

##DotPlot
features <- c("Pdgfra", "C1ql1", "Tnr", "Ust", "Cnksr3", "Ptgds", "Mobp", "Mog", "Mag", "Mbp", "Plp1", "Celf4", "Rbfox3", "Slc17a6", "Slc6a1", "Slc6a5", "Gad2", "Chat","Slc1a2", "Aqp4", "Ly86", "Runx1", "Vtn", "Abcc9", "Flt1","Dnah12")
DotPlot(AR.integrated, features = features) 

#Exclude clusters containing doublets or background
AR.integrated.subset <- subset(AR.integrated,  idents = c("9"), invert=TRUE)
p1 <- DimPlot(AR.integrated.subset, reduction = "umap", group.by = "AR")
plot_grid(p1)　

#Save
saveRDS(AR.integrated.subset,"~/Downloads/R/230908 A36913 merge DEG_before_annotation_after_subset.rds") 

#Label clusters 
new.cluster.ids <- c("Oligodendrocytes","Oligodendrocytes","Oligodendrocytes","Oligodendrocytes",
                     "Oligodendrocytes","Inhibitory neurons","Inhibitory neurons", "Astrocytes",
                     "OPC", "Microglia","Inhibitory neurons","Oligodendrocytes", 
                     "Astrocytes","Excitatory neurons","Excitatory neurons","Oligodendrocytes", 
                     "Excitatory neurons","Excitatory neurons","Excitatory neurons","OL progenitors",
                     "Inhibitory neurons","Excitatory neurons","Astrocytes","Pericytes",
                     "Inhibitory neurons","Endothelial cells","Inhibitory neurons","Ependymal cells")
names(new.cluster.ids) <- levels(AR.integrated.subset)
AR.integrated.subset <- RenameIdents(AR.integrated.subset, new.cluster.ids)

DimPlot(AR.integrated.subset , reduction = "umap", label = FALSE, split.by="AR")　
DimPlot(AR.integrated.subset , reduction = "umap", label = FALSE, 
        cols = c("Oligodendrocytes"="#F8766D","Inhibitory neurons" ="#CD9600", "Astrocytes"="#7CAE00", "Excitatory neurons"="#0CB702", "OPC"="#ED68ED", "OL progenitors"="#00A9FF", "Microglia"="#00BFC4", "Pericytes"="#8494FF","Ependymal cells"="#FF68A1", "Endothelial cells"="#C77CFF" ),
        order = c("Ependymal cells","Endothelial cells","Pericytes","Microglia","OL progenitors","OPC","Astrocytes","Excitatory neurons", "Inhibitory neurons","Oligodendrocyte"))

#DEGs of OLs clusters
AR.integrated.subset$cluster_sample <- paste(AR.integrated.subset$seurat_clusters, AR.integrated.subset$orig.ident, sep = "_")
Idents(AR.integrated.subset) <- "cluster_sample"
DEG_oligo.AD_vs_other <- FindMarkers(AR.integrated.subset, ident.1 = c("0_AD", "1_AD", "2_AD","3_AD", "4_AD", "12_AD", "16_AD" ), ident.2 = c("0_AA", "1_AA", "2_AA","3_AA", "4_AA", "12_AA", "16_AA","0_AB", "1_AB", "2_AB","3_AB", "4_AB", "12_AB", "16_AB","0_AC", "1_AC", "2_AC","3_AC", "4_AC", "12_AC", "16_AC"), logfc.threshold = 0.1)

#GO enrichment analysis of 100 genes elevated in AR-97Q mice
DEG_oligo.AD_vs_other_0.05 <- DEG_oligo.AD_vs_other[DEG_oligo.AD_vs_other$p_val_adj < 0.05
                                                    & DEG_oligo.AD_vs_other$avg_log2FC>
                                                      0.2745,]
DEG.gene_SYMBOLs <- rownames(DEG_oligo.AD_vs_other_0.05)
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

#GO enrichment analysis of 94 genes elevated in AR-97Q mice
DEG_oligo.AD_vs_other_0.05 <- DEG_oligo.AD_vs_other[DEG_oligo.AD_vs_other$p_val_adj < 0.05
                                                    & DEG_oligo.AD_vs_other$avg_log2FC< -0.1
                                                    ,]

DEG.gene_SYMBOLs <- rownames(DEG_oligo.AD_vs_other_0.05)
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
AR.integrated.subset <- readRDS("~/Downloads/R/230908 A36913 merge DEG_before_annotation_after_subset.rds")
DEG_oligo.AD_vs_other <- FindMarkers(AR.integrated.subset, ident.1 = c("0_AD", "1_AD", "2_AD","3_AD", "4_AD", "12_AD", "16_AD" ), ident.2 = c("0_AA", "1_AA", "2_AA","3_AA", "4_AA", "12_AA", "16_AA","0_AB", "1_AB", "2_AB","3_AB", "4_AB", "12_AB", "16_AB","0_AC", "1_AC", "2_AC","3_AC", "4_AC", "12_AC", "16_AC"   ), logfc.threshold = 0.1)

#Annotate the Ensembl gene IDs to gene symbols:
ens <- rownames(DEG_oligo.AD_vs_other)
symbols <- mapIds(org.Mm.eg.db, keys = ens,
                  column = c('SYMBOL'), keytype = 'SYMBOL')
symbols <- symbols[!is.na(symbols)]
symbols <- symbols[match(rownames(DEG_oligo.AD_vs_other), names(symbols))]
rownames(DEG_oligo.AD_vs_other) <- symbols
keep <- !is.na(rownames(DEG_oligo.AD_vs_other))
DEG_oligo.AD_vs_other  <- DEG_oligo.AD_vs_other [keep,]

EnhancedVolcano(DEG_oligo.AD_vs_other, 
                lab = rownames(DEG_oligo.AD_vs_other),
                x ="avg_log2FC", 
                y ="p_val_adj",
                pCutoff = 10e-2, 
                FCcutoff = 0.2,
                pointSize = 1.0,selectLab = c( "9330111N05Rik","Garnl3",'A330015K06Rik', "Kcnma1","Rhoj", "Dpp10","Sgcz","Nrg3", "Cntnap2", "P4ha1" ),
                labSize = 5.0,
                drawConnectors=TRUE, widthConnectors=0.5)  
##########
##########
##########
##########
#Load libraries for Monocle
library(monocle3)
library(SeuratWrappers)
library(cowplot)
library(dplyr)

AR.integrated.OL <- subset(AR.integrated,  idents = c("0", "1","2", "3","4", "8", "12", "16", "20"))

Idents(AR.integrated.OL) <- "seurat_clusters"
DimPlot(AR.integrated.OL)
DefaultAssay(AR.integrated.OL) <- "RNA"
AR.integrated.OL.cds <- as.cell_data_set(AR.integrated.OL)

AR.integrated.OL.cds <- estimate_size_factors(AR.integrated.OL.cds)
rowData(AR.integrated.OL.cds)$gene_short_name <- row.names(rowData(AR.integrated.OL.cds))

AR.integrated.OL.cds <- cluster_cells(cds = AR.integrated.OL.cds, reduction_method = "UMAP", resolution = 0.001)
AR.integrated.OL.cds@clusters@listData[["UMAP"]][["clusters"]]
AR.integrated.OL.cds <- learn_graph(AR.integrated.OL.cds, use_partition = FALSE)
plot_cells(AR.integrated.OL.cds, cell_size = 1)

AR.integrated.OL.cds <- order_cells(AR.integrated.OL.cds)

#Plot pseudotime
plot_cells(AR.integrated.OL.cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)　　

plot_cells(AR.integrated.OL.cds,
           color_cells_by = "AR",
           #group_cells_by = "partition",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5) + facet_wrap(~AR, nrow=2)  

#Gene expression changes in pseudotime
target_genes<-c("Pdgfra")
target_genes<-c("Sox6")
target_genes<-c("Mog")
target_genes<-c("Tnr")
target_genes<-c("Apc")
target_genes<-c("Asic2")
target_genes<-c("Fam155a")

synapse.cds<-AR.integrated.OL.cds[rowData(AR.integrated.OL.cds)$gene_short_name %in% target_genes,
                                  colData(AR.integrated.OL.cds)$seurat_clusters %in% c("0", "1","2", "3","4", "8", "12", "16", "20")]
plot_genes_in_pseudotime(synapse.cds, color_cells_by = "seurat_clusters", min_expr=0.5)


