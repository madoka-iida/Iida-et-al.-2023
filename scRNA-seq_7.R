#The scRNA-seq dataset was processed and analyzed with Seurat version v4.3.1.
#Load libraries
library(Seurat)
library(CellChat)
library(patchwork)
library(NMF)
library(ggalluvial)
library(ComplexHeatmap)

#Load the dataset
seurat_object <- readRDS("WA6.combined_anotated_231120.rds")
Idents(seurat_object) <- "WA6"
A6 <- subset(seurat_object, idents="A6")
DefaultAssay(A6) <- "RNA" 
A6 <- NormalizeData(A6)
cellchat <- createCellChat(object = A6, group.by = "annotated")　

#Load the database of ligand-receptor interaction
CellChatDB <- CellChatDB.mouse 
showDatabaseCategory(CellChatDB)
cellchat@DB <- CellChatDB

#Preprocessing of gene expression data for cell-cell communication analysis 
cellchat <- subsetData(cellchat) 
future::plan("multisession", workers = 4) 

cellchat <- identifyOverExpressedGenes(cellchat, thresh.p = 0.1)　
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat) 

#Communication involving less than 10 cells is filtered.
cellchat <- filterCommunication(cellchat, min.cells = 10)

#Estimation of intercellular signaling at the level of signaling pathways
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents)) 

par(mfrow = c(1,2), xpd=TRUE) 
#Plot the number of intercellular communications and weight/strength
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Interaction weights/strength")

#Visualize intercellular communication for each cluster
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE) 
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}  

#Visualize each signaling pathway
cellchat@netP$pathways 
pathways.show <- c("NCAM")  
levels(cellchat@idents) 
vertex.receiver = c(2, 3, 4) 
netVisual_aggregate(cellchat, signaling = pathways.show, layout="hierarchy", vertex.receiver = vertex.receiver) #plot

#Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

#Outgoing patterns
selectK(cellchat, pattern = "outgoing") 
nPatterns = 4  
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)

#River plot
netAnalysis_river(cellchat, pattern = "outgoing", cutoff=0.4)

#Incoming patterns
selectK(cellchat, pattern = "incoming")
nPatterns = 2 
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)

#River plot
netAnalysis_river(cellchat, pattern = "incoming", cutoff=0.4) ##defaultではcutoff=0.5. 

#Save
saveRDS(cellchat, "cellchat_A6_231121.rds")

###########
#The other sample is analyzed in the same way.
###########

#Load CellChat object of each dataset and merge together
cellchat.W6 <- readRDS("cellchat_W6_231121.rds")
cellchat.A6 <- readRDS("cellchat_A6_231121.rds")
cellchat.W6 <- updateCellChat(cellchat.W6)
cellchat.A6 <- updateCellChat(cellchat.A6)
object.list <- list(W6 = cellchat.W6, A6 = cellchat.A6)
cellchat <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)

#Compare interactions
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2　　

#Differential number of interactions or interaction strength among different cell populations
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, label.edge = TRUE) 
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", label.edge = TRUE) 

#Heatmap based on a merged object
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2　　

#Identify signaling changes associated with one cell group
pdf("signalChnalges_scatterWA6_231121.pdf", width= 15, height = 4) 
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "OL")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "OLpro")
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "OPC")
patchwork::wrap_plots(plots = list(gg1,gg2,gg3)) 
dev.off()

i = 1
#Compare outgoing (or incoming) signaling associated with each cell population 
#Outgoing
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))  

#Incoming
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

#All
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


#Identify the upgulated and down-regulated signaling ligand-receptor pairs
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(5:11),  comparison = c(1, 2), angle.x = 45)
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), angle.x = 45)
netVisual_bubble(cellchat, sources.use = 7, targets.use = c(5:11),  comparison = c(1, 2), angle.x = 45)

gg1 <- netVisual_bubble(cellchat, sources.use = 1, targets.use = c(3,6),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in A6", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 1, targets.use = c(3,6),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in A6", angle.x = 45, remove.isolate = T)
gg1 + gg2

#Save
saveRDS(cellchat, "cellchat_mergeWA6_231121.rds")



