############################################## Step1: QC ##############################################
setwd("~/Documents/NJIT/Rutgers/scRNAseq/")
library(Seurat)
library(DoubletFinder)

nCount_cutoff <- 90000
nFeature_cutoff <- 1000
pMito_cutoff <- 10

hom <- read.table(file = "HOM_MouseHumanSequence.txt", header = T, sep = "\t")
convert_Human_to_Mouse_symbol <- function(x){
  key <- hom$DB.Class.Key[which(hom$Symbol == x)]
  value <- hom$Symbol[which(hom$DB.Class.Key == key & hom$Common.Organism.Name == "mouse, laboratory")]
  return(value)
}

m.s.genes <- sapply(cc.genes.updated.2019$s.genes, convert_Human_to_Mouse_symbol, simplify = T)
m.g2m.genes <- sapply(cc.genes.updated.2019$g2m.genes, convert_Human_to_Mouse_symbol, simplify = T)

# GF
GF <- Read10X_h5(filename = "filtered_feature_bc_matrix_GF.h5")
GF <- CreateSeuratObject(counts = GF)
GF <- RenameCells(object = GF, add.cell.id = "GF")
GF$Condition = "GF"

GF <- PercentageFeatureSet(GF, pattern = "^mt-", col.name = "percent_mito")

GF$QC <- "PASS"
GF$QC[which(GF$nCount_RNA > nCount_cutoff)] <- "High_nCount"
GF$QC[which(GF$nFeature_RNA < nFeature_cutoff)] <- "Low_nFeature"
GF$QC[which(GF$percent_mito > pMito_cutoff)] <- "High_pMito"
GF$QC[which(GF$nCount_RNA > nCount_cutoff & GF$nFeature_RNA < nFeature_cutoff)] <- "High_nCount; Low_nFeature"
GF$QC[which(GF$nCount_RNA > nCount_cutoff & GF$percent_mito > pMito_cutoff)] <- "High_nCount; High_pMito"
GF$QC[which(GF$nFeature_RNA < nFeature_cutoff & GF$percent_mito > pMito_cutoff)] <- "Low_nFeature; High_pMito"
GF$QC[which(GF$nCount_RNA > nCount_cutoff & GF$nFeature_RNA < nFeature_cutoff & GF$percent_mito > pMito_cutoff)] <- "High_nCount; Low_nFeature; High_pMito"

temp <- GF
temp <- NormalizeData(temp, verbose = FALSE)
temp <- FindVariableFeatures(temp, verbose = FALSE)
temp <- ScaleData(temp, verbose = FALSE)
temp <- RunPCA(temp, npcs = 10, verbose = FALSE)
temp <- RunUMAP(temp, dims = 1:10, verbose = FALSE)

temp <- CellCycleScoring(temp, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)
GF$S.Score <- temp$S.Score
GF$G2M.Score <- temp$G2M.Score
GF$Phase <- temp$Phase

sweep.res <- paramSweep(temp, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
best_pK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
nExp_poi <- round(0.05 * nrow(temp@meta.data))
temp <- doubletFinder(temp, PCs = 1:10, pN = 0.25, pK = best_pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
GF$"DF_classifications" <- temp[[paste("DF.classifications_0.25", best_pK, nExp_poi, sep = "_")]]

Immune_cells <- CellSelector(FeaturePlot(temp, features = "Ptprc"))
Epithelial_cells  <- CellSelector(FeaturePlot(temp, features = "Muc3"))

GF$Annotation <- "Paneth cells"
GF$Annotation[which(colnames(GF) %in% Immune_cells)] <- "Immune cells"
GF$Annotation[which(colnames(GF) %in% Epithelial_cells)] <- "Intestinal epithelial cells"

saveRDS(object = GF, file = "GF.rds")

# SPF
SPF <- Read10X_h5(filename = "filtered_feature_bc_matrix_BL6M.h5")
SPF <- CreateSeuratObject(counts = SPF)
SPF <- RenameCells(object = SPF, add.cell.id = "SPF")
SPF$Condition = "SPF"

SPF <- PercentageFeatureSet(SPF, pattern = "^mt-", col.name = "percent_mito")

SPF$QC <- "PASS"
SPF$QC[which(SPF$nCount_RNA > nCount_cutoff)] <- "High_nCount"
SPF$QC[which(SPF$nFeature_RNA < nFeature_cutoff)] <- "Low_nFeature"
SPF$QC[which(SPF$percent_mito > pMito_cutoff)] <- "High_pMito"
SPF$QC[which(SPF$nCount_RNA > nCount_cutoff & SPF$nFeature_RNA < nFeature_cutoff)] <- "High_nCount; Low_nFeature"
SPF$QC[which(SPF$nCount_RNA > nCount_cutoff & SPF$percent_mito > pMito_cutoff)] <- "High_nCount; High_pMito"
SPF$QC[which(SPF$nFeature_RNA < nFeature_cutoff & SPF$percent_mito > pMito_cutoff)] <- "Low_nFeature; High_pMito"
SPF$QC[which(SPF$nCount_RNA > nCount_cutoff & SPF$nFeature_RNA < nFeature_cutoff & SPF$percent_mito > pMito_cutoff)] <- "High_nCount; Low_nFeature; High_pMito"

temp <- SPF
temp <- NormalizeData(temp, verbose = FALSE)
temp <- FindVariableFeatures(temp, verbose = FALSE)
temp <- ScaleData(temp, verbose = FALSE)
temp <- RunPCA(temp, npcs = 10, verbose = FALSE)
temp <- RunUMAP(temp, dims = 1:10, verbose = FALSE)

temp <- CellCycleScoring(temp, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)
SPF$S.Score <- temp$S.Score
SPF$G2M.Score <- temp$G2M.Score
SPF$Phase <- temp$Phase

sweep.res <- paramSweep(temp, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
best_pK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
nExp_poi <- round(0.05 * nrow(temp@meta.data))
temp <- doubletFinder(temp, PCs = 1:10, pN = 0.25, pK = best_pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
SPF$"DF_classifications" <- temp[[paste("DF.classifications_0.25", best_pK, nExp_poi, sep = "_")]]

Immune_cells <- CellSelector(FeaturePlot(temp, features = "Ptprc"))
Epithelial_cells  <- CellSelector(FeaturePlot(temp, features = "Muc3"))

SPF$Annotation <- "Paneth cells"
SPF$Annotation[which(colnames(SPF) %in% Immune_cells)] <- "Immune cells"
SPF$Annotation[which(colnames(SPF) %in% Epithelial_cells)] <- "Intestinal epithelial cells"

saveRDS(object = SPF, file = "SPF.rds")

########################################## Step2: Integration #########################################
setwd("~/Documents/NJIT/Rutgers/scRNAseq/")
library(Seurat)
options(future.globals.maxSize = 24000 * 1024^2)

GF <- readRDS(file = "GF.rds")
SPF <- readRDS(file = "SPF.rds")

GF <- subset(GF, QC == "PASS" & DF_classifications == "Singlet" & Annotation == "Paneth cells")
SPF <- subset(SPF, QC == "PASS" & DF_classifications == "Singlet" & Annotation == "Paneth cells")

data <- merge(GF, SPF)
data <- JoinLayers(data, assay = "RNA")

data <- SCTransform(object = data, vst.flavor = "v2", do.correct.umi = T, vars.to.regress = c("S.Score", "G2M.Score"))
data <- RunPCA(data)
ElbowPlot(data)
data <- RunUMAP(data, reduction = "pca", dims = 1:5)

DimPlot(data, reduction = "umap", group.by = "Condition")
DimPlot(data, reduction = "umap", group.by = "Phase")
FeaturePlot(object = data, reduction = "umap", features = "Mki67", pt.size = 0.5)

saveRDS(object = data, file = "Integrated_Data.rds")

########################################## Step3: Differential ######################################## 
setwd("~/Documents/NJIT/Rutgers/scRNAseq/")
library(Seurat)
library(EnhancedVolcano)
library(dplyr)
library(fgsea)
library(AUCell)

data <- readRDS(file = "Integrated_Data.rds")

tiff(filename = "~/Documents/NJIT/Rutgers/scRNAseq/Fgfbp1.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(object = data, reduction = "umap", features = "Fgfbp1", pt.size = 0.3, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "~/Documents/NJIT/Rutgers/scRNAseq/Mki67.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(object = data, reduction = "umap", features = "Mki67", pt.size = 0.3, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "~/Documents/NJIT/Rutgers/scRNAseq/Lyz1.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(object = data, reduction = "umap", features = "Lyz1", pt.size = 0.3, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "~/Documents/NJIT/Rutgers/scRNAseq/Mptx2.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(object = data, reduction = "umap", features = "Mptx2", pt.size = 0.3, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "~/Documents/NJIT/Rutgers/scRNAseq/Cdc42.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(object = data, reduction = "umap", features = "Cdc42", pt.size = 0.3, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "~/Documents/NJIT/Rutgers/scRNAseq/Defa27.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(object = data, reduction = "umap", features = "Defa27", pt.size = 0.3, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "~/Documents/NJIT/Rutgers/scRNAseq/Cdc42.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(object = data, reduction = "umap", features = "Cdc42", pt.size = 0.3, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "~/Documents/NJIT/Rutgers/scRNAseq/Rac1.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(object = data, reduction = "umap", features = "Rac1", pt.size = 0.3, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "~/Documents/NJIT/Rutgers/scRNAseq/Rhoa.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(object = data, reduction = "umap", features = "Rhoa", pt.size = 0.3, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "~/Documents/NJIT/Rutgers/scRNAseq/Olfm4.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(object = data, reduction = "umap", features = "Olfm4", pt.size = 0.3, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()


tiff(filename = "~/Documents/NJIT/Rutgers/scRNAseq/Cdc42 (original color).tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(object = data, reduction = "umap", features = "Cdc42", pt.size = 0.3)
dev.off()

tiff(filename = "~/Documents/NJIT/Rutgers/scRNAseq/Rac1 (original color).tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(object = data, reduction = "umap", features = "Rac1", pt.size = 0.3)
dev.off()

tiff(filename = "~/Documents/NJIT/Rutgers/scRNAseq/Rhoa (original color).tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(object = data, reduction = "umap", features = "Rhoa", pt.size = 0.3)
dev.off()

tiff(filename = "~/Documents/NJIT/Rutgers/scRNAseq/Olfm4 (original color).tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(object = data, reduction = "umap", features = "Olfm4", pt.size = 0.3)
dev.off()

library(msigdbr)
counts <- GetAssayData(object = data, assay = "SCT", layer = "data")
genesets <- msigdbr(species = "Mus musculus", collection = "H")
gene_sets <- split(x = genesets$gene_symbol, f = genesets$gs_name)
pathway <- gene_sets$HALLMARK_OXIDATIVE_PHOSPHORYLATION
data$HALLMARK_OXIDATIVE_PHOSPHORYLATION <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = pathway)))
pathway <- gene_sets$HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY
data$HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = pathway)))

tiff(filename = "~/Documents/NJIT/Rutgers/scRNAseq/HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(data, reduction = "umap", features = "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", , cols = c("lightgrey", "#FFE680", "#FF6B6B")) + theme(plot.title = element_text(size = 10))
dev.off()

data <- FindNeighbors(data, reduction = "pca")
data <- FindClusters(data, resolution = 0.2)
DimPlot(data, reduction = "umap", group.by = "seurat_clusters")

data$group <- "Progenitor Paneth cell"
data$group[which(data$seurat_clusters == 1)] <- "GF Paneth cells"
data$group[which(data$seurat_clusters %in% c(3,4))] <- "SPF Paneth cells"
data$group <- factor(x = data$group, levels = c("Progenitor Paneth cell", "SPF Paneth cells", "GF Paneth cells"))

res <- FindMarkers(object = data, ident.1 = "SPF", ident.2 = "GF", group.by = "group")

EnhancedVolcano(toptable = res, 
                lab = rownames(res), 
                x = "avg_log2FC",
                y = "p_val",
                title = "Volcano plot",
                subtitle = "SPF vs. GF",
                pCutoff = 0.05, 
                pCutoffCol = "p_val_adj", 
                FCcutoff = 2, 
                cutoffLineType = "dashed")

########################################### Step4: Potency ###########################################
setwd("~/Documents/NJIT/Rutgers/scRNAseq/")
library(Seurat)
library(CytoTRACE2)

data <- readRDS(file = "Integrated_Data.rds")
data <- cytotrace2(data, species = "mouse", is_seurat = T, slot_type = "count", ncores = 1)

df <- data@meta.data
ggplot(data = df, mapping = aes(x = Condition, y = CytoTRACE2_Score, fill = Condition)) +
  geom_violin(width = 1) +
  geom_boxplot(width = 0.1) +
  theme_classic()

t.test(x = df$CytoTRACE2_Score[which(df$Condition == "GF")],
       y = df$CytoTRACE2_Score[which(df$Condition == "SPF")])

df2 <- as.data.frame(table(df$Condition, df$CytoTRACE2_Potency))
colnames(df2) <- c("Condition", "Potency", "Counts")
ggplot(data = df2, mapping = aes(x = Condition, y = Counts, fill = Potency)) +
  geom_bar(position="fill", stat="identity") +
  xlab(label = "Condition") +
  ylab(label = "Proportion") +
  theme_classic()

DimPlot(data, reduction = "umap", group.by = "CytoTRACE2_Potency")
FeaturePlot(data, features = "Gpx2")

# subset undifferentiated cells
data2 <- subset(data, CytoTRACE2_Potency %in% c("Unipotent", "Oligopotent", "Multipotent"))

cdc42_downstream_mouse <- c("Cdc42",
                            "Pak1", "Pak2", "Wasl", "Wasf2",
                            "Arpc1b", "Arpc3", "Iqgap1", "Rac1")

tiff(filename = "~/Documents/NJIT/Rutgers/scRNAseq/Dotplot.tiff", width = 900, height = 320, res = 100)
DotPlot(data2, assay = "SCT", features = cdc42_downstream_mouse, group.by = "Condition", col.min = 0) +
  RotatedAxis() +
  labs(title = "Cdc42 Downstream Genes by Condition")
dev.off()

# correlation
expr_matrix <- GetAssayData(data, assay = "SCT", slot = "data")
cors <- apply(expr_matrix, 1, function(gene_expr) {
  cor(gene_expr, data$CytoTRACE2_Score, method = "pearson")
})

cors["Cdc42"]
cors["Stmn1"]

WriteXLS::WriteXLS(x = as.data.frame(cors), ExcelFileName = "Pearson correlation between Cdc42 and CytoTRACE2_Score.xlsx", row.names = T)

df <- data.frame(Cdc42 = expr_matrix["Cdc42",],
                 Potency = data$CytoTRACE2_Score)

ggplot(data = df, mapping = aes(x = Cdc42, y = Potency)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_continuous(transform = "log1p") +
  scale_y_continuous(transform = "log1p") +
  theme_classic()

# heatmap
library(pheatmap)

gene_list <- c("Stmn1","Olfm4","Fgfbp1","Cdc42","Mki67","Rhoa","Lgr5","Rac1","Keap1","Dusp3","Mptx2","Defa27","Mmp7","Lyz1")
pathway_list <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", 
                  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY")

draw_heatmap <- function(seurat_obj, cells_to_use, title_text) {
  cytotrace_scores <- seurat_obj$CytoTRACE2_Score[cells_to_use]
  ordered_cells <- names(sort(cytotrace_scores, decreasing = TRUE))
  
  expr_matrix <- GetAssayData(seurat_obj, layer = "data")
  expr_genes <- expr_matrix[gene_list, ordered_cells]
  
  pathway_scores <- t(seurat_obj@meta.data[ordered_cells, pathway_list])
  rownames(pathway_scores) <- pathway_list
  
  combined_matrix <- rbind(as.matrix(expr_genes), as.matrix(pathway_scores))
  
  annotation_col <- data.frame(CytoTRACE = cytotrace_scores[ordered_cells])
  rownames(annotation_col) <- ordered_cells
  
  annotation_colors <- list(CytoTRACE = colorRampPalette(c("white", "red"))(100))
  
  pheatmap(combined_matrix,
           scale = "row",
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           show_rownames = TRUE,
           show_colnames = FALSE,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           annotation_col = annotation_col,
           annotation_colors = annotation_colors,
           main = title_text, 
           gaps_row = length(gene_list))
}

tiff(filename = "~/Documents/NJIT/Rutgers/scRNAseq/Heatmap.tiff", width = 2000, height = 1260, res = 200)
draw_heatmap(data, colnames(data), "All Cells")
dev.off()

draw_heatmap(data, WhichCells(data, expression = Condition == "GF"), "GF Cells")
draw_heatmap(data, WhichCells(data, expression = Condition == "SPF"), "SPF Cells")

#### Step5: Trajectory ####
setwd("~/Documents/NJIT/Rutgers/scRNAseq/")
library(Seurat)
library(monocle3)
library(ggplot2)

data <- readRDS(file = "Integrated_Data.rds")

expr_matrix <- GetAssayData(data, slot = "data")
cell_metadata <- data@meta.data
gene_annotation <- as.data.frame(rownames(data), row.names = rownames(data))
colnames(gene_annotation) <- "gene_short_name"

cds <- new_cell_data_set(expr_matrix, 
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 10)
cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA")

reducedDims(cds)$UMAP <- data@reductions$umap@cell.embeddings
cds <- cluster_cells(cds, resolution = 0.0015)

cds <- learn_graph(cds)
cds <- order_cells(cds)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           cell_size = 0.5, 
           trajectory_graph_color = "grey30")

fit_res <- fit_models(cds, model_formula_str = "~ pseudotime")
coef_res <- coefficient_table(fit_res)
coef_res <- subset(coef_res, term == "pseudotime")
coef_res <- coef_res[,c("gene_id", "num_cells_expressed", "status", "term", "estimate", "std_err", "test_val", "p_value", "normalized_effect", "q_value")]

WriteXLS::WriteXLS(x = coef_res, ExcelFileName = "coefficient table.xlsx")

condition_genes <- subset(coef_res, term == "ConditionSPF" & q_value < 0.05 & abs(estimate) > 2)
plot_genes_in_pseudotime(cds[c("Reg3g","Defa35","Defa27"),], color_cells_by = "Condition")
plot_genes_in_pseudotime(cds[c("Defa23","Defa36"),], color_cells_by = "Condition")
plot_genes_in_pseudotime(cds["Cdc42",], color_cells_by = "Condition")
plot_genes_in_pseudotime(cds["Olfm4",], color_cells_by = "Condition")
plot_genes_in_pseudotime(cds["Cdc42",], color_cells_by = "Condition")

cds_gf <- cds[, colData(cds)$Condition == "GF"]
cds_spf <- cds[, colData(cds)$Condition == "SPF"]

p2 <- plot_genes_in_pseudotime(cds_gf["Mptx2", ], color_cells_by = "pseudotime") + 
  ylim(c(-1,1600)) + 
  ggtitle("GF")
p1 <- plot_genes_in_pseudotime(cds_spf["Mptx2", ], color_cells_by = "pseudotime") + 
  ylim(c(-1,1600)) +  
  ggtitle("SPF")

tiff(filename = "~/Documents/NJIT/Rutgers/scRNAseq/Mptx2 trajectory.tiff", width = 1200, height = 240, res = 100)
plot_grid(p1, p2, ncol = 2, label_size = 12)
dev.off()

plot_genes_in_pseudotime(cds_gf["Defa27", ], color_cells_by = "pseudotime") + 
  ylim(c(0,250)) + 
  ggtitle("GF")
plot_genes_in_pseudotime(cds_spf["Defa27", ], color_cells_by = "pseudotime") + 
  ylim(c(0,250)) +  
  ggtitle("SPF")

plot_genes_in_pseudotime(cds_gf["Defa35", ], color_cells_by = "pseudotime") + 
  ylim(c(0,750)) + 
  ggtitle("GF")
plot_genes_in_pseudotime(cds_spf["Defa35", ], color_cells_by = "pseudotime") + 
  ylim(c(0,750)) +  
  ggtitle("SPF")

plot_genes_in_pseudotime(cds_gf["Reg3g", ], color_cells_by = "pseudotime") + 
  ylim(c(0,400)) + 
  ggtitle("GF")
plot_genes_in_pseudotime(cds_spf["Reg3g", ], color_cells_by = "pseudotime") + 
  ylim(c(0,400)) +  
  ggtitle("SPF")

plot_genes_in_pseudotime(cds_gf["Defa34", ], color_cells_by = "pseudotime") + 
  ylim(c(0,250)) + 
  ggtitle("GF")
plot_genes_in_pseudotime(cds_spf["Defa34", ], color_cells_by = "pseudotime") + 
  ylim(c(0,250)) +  
  ggtitle("SPF")

plot_genes_in_pseudotime(cds_gf["Mptx2", ], color_cells_by = "pseudotime") + 
  ylim(c(0,1500)) + 
  ggtitle("GF")
plot_genes_in_pseudotime(cds_spf["Mptx2", ], color_cells_by = "pseudotime") + 
  ylim(c(0,1500)) +  
  ggtitle("SPF")

