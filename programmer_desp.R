library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(tximport)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)

dir <- ("/projectnb/bf528/project_4_scrnaseq/GSM2230760__salmon_quant")

files <- file.path(dir, "alevin", "quants_mat.gz")
file.exists(files)
txi <- tximport(files, type="alevin")

# need to remove the number after period ENSXXXX.number
ensmbl_names <- txi$counts@Dimnames[[1]]
ensmbl.clean <- sapply(strsplit(as.character(ensmbl_names), "\\."), "[[", 1)

#keytypes(org.Hs.eg.db)
cheat.sheet <- AnnotationDbi::select(org.Hs.eg.db, ensmbl.clean, 
                                     "GENENAME", "ENSEMBL")

length(cheat.sheet$ENSEMBL) #60465
length(ensmbl_names) #60233
# this is very annoying there are duplicates in the cheat.sheet
# need to find duplicates and remove one

cheat.sheet <- cheat.sheet[!duplicated(cheat.sheet$ENSEMBL), ]
length(cheat.sheet$ENSEMBL) #60233
# ok maybe I was being dramatic it wasn't that big a deal

# not all ensembl.ids have symbols though
# if GENENAME is.na use ensembl id?

##ifelse(expression, true, false)
new.rownames <- ifelse(is.na(cheat.sheet$GENENAME), cheat.sheet$ENSEMBL, 
                       cheat.sheet$GENENAME)

length(new.rownames) #60233

txi$counts@Dimnames[[1]] <- new.rownames

ctrl <- CreateSeuratObject(counts = txi$counts, min.features = 100, project=".")
dim(ctrl)
#head(ctrl@assays$RNA@counts[1:3, 1:3])
View(ctrl@meta.data)
#nCount_RNA number of UMIs per cell
#nFeature_RNA number of genes detected per cell

# Add number of genes per UMI for each cell to metadata
ctrl$log10GenesPerUMI <- log10(ctrl$nFeature_RNA) / log10(ctrl$nCount_RNA)

# Mitochondrial QC metrics
ctrl[["percent.mt"]] <- PercentageFeatureSet(ctrl, pattern = "^MT-")
#View(ctrl@meta.data)

# Visualise QC metrics as violin plot
jpeg("violin_plot.jpeg")
VlnPlot(ctrl, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
dev.off()

plot1 <- FeatureScatter(ctrl, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ctrl, feature1 = "nCount_RNA", feature2 ="nFeature_RNA")
#plot3 <- FeatureScatter(ctrl,feature1="log10GenesPerUMI",feature2="percent.mt")

oldpar<-par(no.readonly=TRUE)

jpeg("feature_relationships1.jpeg")
plot1
dev.off()

jpeg("feature_relationships2.jpeg")
plot2
dev.off()

ctrl <- subset(ctrl, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 &
                 percent.mt <5) #standard threshold is mt 5%

ctrl <- NormalizeData(ctrl, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
#saved in ctrl[["RNA"]]@data

ctrl <- FindVariableFeatures(ctrl)

var.feat.plt <- VariableFeaturePlot(ctrl)
#from here we see that about 50 are the most variable

# Identify the 50 most highly variable genes
top50 <- head(VariableFeatures(ctrl), 50)
top20 <- head(VariableFeatures(ctrl), 20)
top15 <- head(VariableFeatures(ctrl), 15)
top10 <- head(VariableFeatures(ctrl), 10)
top5 <- head(VariableFeatures(ctrl), 5)
as.list(top15)
write.csv(top15, "top15.csv")

jpeg("variable_features.jpeg")
var.feat.plt
dev.off()

# scaling data
all.genes <- rownames(ctrl)
ctrl <- ScaleData(ctrl, features = all.genes)
# adding covariates, this one took longer to run
ctrl2 <- ScaleData(ctrl, features = all.genes, vars.to.regress = 
                     c("percent.mt","nFeature_RNA"), assay = "RNA")

# linear dimensional reduction
ctrl <- RunPCA(ctrl, features = VariableFeatures(object = ctrl))
ctrl2 <- RunPCA(ctrl2, features = VariableFeatures(object = ctrl))

#Visualize top genes associated with reduction components
jpeg("dimensional_redution_genes.jpeg")
VizDimLoadings(ctrl, dims = 1:2, reduction = "pca") 
dev.off()


dim.plt1 <- DimPlot(ctrl, reduction = "pca")
dim.plt2 <- DimPlot(ctrl2, reduction = "pca")
jpeg("dim_plots.jpeg")
dim.plt1
dev.off()

#print(ctrl[["pca"]], dims = 1:5, nfeatures = 5)
# which are significant though?
# didn't do the heatmaps bc I don't think they're helpful for anything

ctrl <- JackStraw(ctrl, num.replicate = 100)
ctrl <- ScoreJackStraw(ctrl, dims = 1:20)
jpeg("JackStraw_nocovar_20.jpeg")
JackStrawPlot(ctrl, dims = 1:20)
dev.off()
jpeg("JackStraw_nocovar_15.jpeg")
JackStrawPlot(ctrl, dims = 1:15)
dev.off()
# amount of variance explained by each PC
jpeg("elbow_nocovar.jpeg")
ElbowPlot(ctrl)
dev.off()

#ctrl2
ctrl2 <- JackStraw(ctrl2, num.replicate = 100)
ctrl2 <- ScoreJackStraw(ctrl2, dims = 1:20)
jpeg("JackStraw_covar_20.jpeg")
JackStrawPlot(ctrl2, dims = 1:20)
dev.off()
jpeg("JackStraw_covar_15.jpeg")
JackStrawPlot(ctrl2, dims = 1:15)
dev.off()

# amount of variance explained by each PC
jpeg("elbow_covar.jpeg")
ElbowPlot(ctrl2)
dev.off()

# will choose 7 PCs
ctrl <- FindNeighbors(ctrl, dims = 1:9)
ctrl2 <- FindNeighbors(ctrl2, dims = 1:6) #if more, I only get 12 clusters

ctrl <- FindClusters(ctrl, resolution = 0.5)
ctrl2 <- FindClusters(ctrl2, resolution = 0.5)

# ctrl and ctrl2 assign barcodes to different clusters?!
head(Idents(ctrl), 10) 
head(Idents(ctrl2), 10)

# tabulate number of cells present in each cluster
counts1 <- table(ctrl$seurat_clusters)
counts2 <- table(ctrl2$seurat_clusters)
prop.table(counts1)
prop.table(counts2)
######################################
##########     analyst      ##########
######################################

# Identify markers for each cluster
#1
cluster4.markers <- FindMarkers(ctrl, ident.1 = 4, min.pct = 0.25)
cluster4.markers <- cluster4.markers[which(cluster4.markers$p_val_adj<0.05),]
clu4.tan <- cluster4.markers %>% top_n(n = 10, wt = avg_log2FC)

#2
cluster2.markers <- FindMarkers(ctrl, ident.1 = 2, min.pct = 0.25)
cluster2.markers <- cluster2.markers[which(cluster2.markers$p_val_adj<0.05),]
clu2.tan <- cluster2.markers %>% top_n(n = 5, wt = avg_log2FC)
############ chicle y le sigo a esto

ctrl.markers <- FindAllMarkers(ctrl, only.pos = TRUE, min.pct = 0.25, 
                               logfc.threshold = 0.25)


# need to add gene symbol
cheater <- AnnotationDbi::select(org.Hs.eg.db, ctrl.markers$gene,
                                 "SYMBOL", "GENENAME")

ctrl.markers <- merge(ctrl.markers, cheater, by.x = "gene", by.y = "GENENAME")
ctrl.markers <- ctrl.markers[which(ctrl.markers$p_val_adj<0.05),]
top.markers <-ctrl.markers %>% group_by(cluster) %>% top_n(n = 5, wt=avg_log2FC)



top.markers <- as.data.frame(top.markers)
write.csv(top.markers, "top_markers5.csv")


#############

jpeg("vlnBETA.jpeg")
VlnPlot(ctrl, features = c("insulin"))
dev.off()

jpeg("featBETA.jpeg")
FeaturePlot(ctrl, features = c("insulin"))
dev.off()

jpeg("vlnDUCTAL.jpeg")
VlnPlot(ctrl, features =  c("keratin 19"))
dev.off()

jpeg("featDUCTAL.jpeg")
FeaturePlot(ctrl, features = c("keratin 19"))
dev.off()

jpeg("vlnACINAR.jpeg")
VlnPlot(ctrl, features = c("carboxypeptidase A1"))
dev.off()

jpeg("featacinarL.jpeg")
FeaturePlot(ctrl, features = c("carboxypeptidase A1"))
dev.off()

jpeg("vlnALPHA.jpeg")
VlnPlot(ctrl, features = c("crystallin beta A2"))
dev.off()

jpeg("featALPHA.jpeg")
FeaturePlot(ctrl, features = c("crystallin beta A2"))
dev.off()

jpeg("vlnDELTA.jpeg")
VlnPlot(ctrl, features = c("proteoglycan 4"))
dev.off()

jpeg("featDELTA.jpeg")
FeaturePlot(ctrl, features = c("proteoglycan 4"))
dev.off()

jpeg("vlnSTELLATE.jpeg")
VlnPlot(ctrl, features = c("platelet derived growth factor receptor beta"))
dev.off()

jpeg("featSTELLATE.jpeg")
FeaturePlot(ctrl, features = c("platelet derived growth factor receptor beta"))
dev.off()

jpeg("vlnDUCTAL2.jpeg")
VlnPlot(ctrl, features = c("C-reactive protein"))
dev.off()

jpeg("featDUCTAL2.jpeg")
FeaturePlot(ctrl, features = c("C-reactive protein"))
dev.off()

jpeg("vlnVASC.jpeg")
VlnPlot(ctrl, features = c("platelet and endothelial cell adhesion molecule 1"))
dev.off()

jpeg("featVASC.jpeg")
FeaturePlot(ctrl, features = c("platelet and endothelial cell adhesion molecule 1"))
dev.off()

jpeg("vln0.jpeg")
VlnPlot(ctrl, features = c("adenylate kinase 3"))
dev.off()

jpeg("feat0.jpeg")
FeaturePlot(ctrl, features = c("adenylate kinase 3"))
dev.off()


top10.markers <- ctrl.markers %>% group_by(cluster) %>% top_n(n = 10, 
                                                              wt = avg_log2FC)
DoHeatmap(ctrl, features = top10.markers$gene) + NoLegend()

jpeg("features.jpeg")
FeaturePlot(ctrl, features = "insulin")
dev.off()
top_genes <- as.vector(top50)
top15_genes <- as.vector(top15)
write.csv(ctrl.markers, "ctrl_markers.csv")
write.csv(top50, "top50_genes.csv")

#########################################
new.cluster.ids <- c("0", "beta", "beta", "ductal", "4", "gamma", "acinar", 
                     "alpha", "8", "delta", "stellate", "ductal", "vascular")
names(new.cluster.ids) <- levels(ctrl)
ctrl <- RenameIdents(ctrl, new.cluster.ids)


#non-linear dimensional reduction UMAP
ctrl <- RunUMAP(ctrl, dims = 1:9)



jpeg("UMAP.jpeg")
DimPlot(ctrl, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

#Non-linear dimensional reduction 
ctrl <- RunTSNE(ctrl, dims.use = 1:9 ,reduction.use = "pca", dim_embed = 2)


jpeg("TSNE.jpeg")
DimPlot(ctrl, reduction="tsne", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

