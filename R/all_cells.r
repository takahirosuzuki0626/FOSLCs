#==============================================================================================#

#intiial settings

#==============================================================================================#
# name pattern
run_id <- "all"

# create directories
dir.create(paste0("./out/", run_id))
dir.create(paste0("./data/", run_id))
out_dir <- paste0("./out/", run_id, "/")
data_dir <- paste0("./data/", run_id, "/")

#load utilities
source("R/util.r")

#==============================================================================================#

# load data

#==============================================================================================#
data_paths <- c(
    E10 = "./data/mapped/E10",
    E11 = "./data/mapped/E11",
    E12 = "./data/mapped/E12",
    E13 = "./data/mapped/E13",
    E14 = "./data/mapped/E14",
    iD6 = "./data/mapped/induced_D6"
)

data_list <- NULL
for(i in seq(length(data_paths))){
    data_list <- c(
        data_list,
        list(scCountProc(data_path = data_paths[i], sample_name = names(data_paths)[i], normalization = FALSE))
    )
}
names(data_list) <- names(data_paths)
#==============================================================================================#

# QC

# using QC.r
#==============================================================================================#
dir.create(paste0(out_dir, "QC"))

qc_plot_list <- lapply(data_list, QC)

for(i in seq(length(data_paths))){
    data_list[[i]][["percent.mt"]] <- qc_plot_list[[i]]$percent.mt

    p <- plot_grid(qc_plot_list[[i]]$violin_plot, qc_plot_list[[i]]$scatter_plot, ncol = 1, nrow = 2)
    ggsave(p, file = paste0(out_dir, "QC/", names(data_list)[i], ".pdf"), width = 21, height = 14)
}

# cut-off is manually decided by violin plot &scatter plot
data_list[["E10"]] <- subset(data_list[["E10"]], subset = nFeature_RNA > 750 & nFeature_RNA < 6000 & nCount_RNA < 40000 & percent.mt < 25)
data_list[["E11"]] <- subset(data_list[["E11"]], subset = nFeature_RNA > 750 & nFeature_RNA < 8000 & nCount_RNA < 60000 & percent.mt < 25)
data_list[["E12"]] <- subset(data_list[["E12"]], subset = nFeature_RNA > 750 & nFeature_RNA < 6500 & nCount_RNA < 50000 & percent.mt < 25)
data_list[["E13"]] <- subset(data_list[["E13"]], subset = nFeature_RNA > 750 & nFeature_RNA < 6500 & nCount_RNA < 50000 & percent.mt < 25)
data_list[["E14"]] <- subset(data_list[["E14"]], subset = nFeature_RNA > 750 & nFeature_RNA < 5500 & nCount_RNA < 50000 & percent.mt < 25)
data_list[["iD6"]] <- subset(data_list[["iD6"]], subset = nFeature_RNA > 750 & nFeature_RNA < 5000 & nCount_RNA < 19000 & percent.mt < 25)

#==============================================================================================#

# Normalization

#==============================================================================================#
for(i in seq(length(data_paths))){
    data_list[[i]] <- NormalizeData(object = data_list[[i]], verbose = FALSE)
    data_list[[i]] <- FindVariableFeatures(object = data_list[[i]], selection.method = "vst", nfeatures = 1000, verbose = FALSE)
}

#==============================================================================================#

# Data Integration

#==============================================================================================#
data.anchors <- FindIntegrationAnchors(object.list = data_list, dims = 1:30)
all_genes <- data.anchors@object.list[[1]]@assays$RNA@counts@Dimnames[[1]]
data.integrated <- IntegrateData(anchorset = data.anchors, dims = 1:30, features.to.integrate = all_genes)

age <- data.integrated@meta.data$sample
age %>%
    replace(which(date == "E10"), 10) %>%
    replace(which(date == "E11"), 11) %>%
    replace(which(date == "E12"), 12) %>%
    replace(which(date == "E13"), 13) %>%
    replace(which(date == "E14"), 14) %>%
    replace(which(date == "iD6"), 12) -> age

vivo.vitro <- data.integrated@meta.data$sample
vivo.vitro %>%
    replace(which(vivo.vitro == "E10"), "vivo") %>%
    replace(which(vivo.vitro == "E11"), "vivo") %>%
    replace(which(vivo.vitro == "E12"), "vivo") %>%
    replace(which(vivo.vitro == "E13"), "vivo") %>%
    replace(which(vivo.vitro == "E14"), "vivo") %>%
    replace(which(vivo.vitro == "iD6"), "vitro") -> vivo.vitro

data.integrated@meta.data <- cbind(data.integrated@meta.data, age, vivo.vitro)

save(data.integrated, file = paste0(data_dir, "data.integrated.QC.RData"))

#==============================================================================================#

# Cell Cycle Regression

#==============================================================================================#
dir.create(paste0(out_dir, "ccReg_results"))

data.integrated <- ccReg(data = data.integrated)
save(data.integrated, file = paste0(out_dir, "ccReg_results/data.integrated.QC.ccreg.RData"))

#==============================================================================================#

# Clustering & UMAP (Figure2A)

#==============================================================================================#
dims <- 15
resolution <- 0.2

data.integrated <- RunPCA(object = data.integrated, npcs = dims, verbose = FALSE)
data.integrated <- RunUMAP(object = data.integrated, reduction = "pca", dims = 1:dims)
data.integrated <- FindNeighbors(object = data.integrated, dims = 1:dims)
data.integrated <- FindClusters(object = data.integrated, algorithm = 1, resolution = resolution)

save(data.integrated, file = paste0(data_dir, "data.integrated.QC.ccreg.clust.RData"))

umap <- DimPlot(object = object, reduction = "umap", split.by = "sample", pt.size = 1)
ggsave(umap, file = paste0(out_dir, "by.sample_umap.pdf"), height = 7, width = 7 * nsample)

#==============================================================================================#

# known marker expression (Figure 2B;Supplementary Figure 6B, D)

#==============================================================================================#
dir.create(paste0(out_dir, "known_marker"))

#Violin Plot
multiVlnPlot(object = data.integrated, features = MARKERS$pan, title = "Pan Markers", out_file = paste0(out_dir, "known_marker/pan_vln.pdf"))

#Feature Plot
pdf(paste0(out_dir, "known_marker/granulosa_featureplot.pdf"), width = 28, height = 14)
FeaturePlot(object = data.integrated, features = MARKERS$granulosa, min.cutoff = 0, max.cutoff = c(4.0, 4.0, 3.0), cols = c("gray", "red"), ncol = 4)
dev.off()

pdf(paste0(out_dir, "known_marker/progenitor_featureplot.pdf"), width = 21, height = 7)
FeaturePlot(object = data.integrated, features = MARKERS$progenitor, min.cutoff = 0, max.cutoff = c(3.5, 2.5, 1.8), cols = c("gray", "red"), ncol = 3)
dev.off()

pdf(paste0(out_dir, "known_marker/stromal_featureplot.pdf"), width = 28, height = 21)
FeaturePlot(object = data.integrated, features = MARKERS$stromal, min.cutoff = 0, max.cutoff = c(2.0, 2.0, 3, 3.0), cols = c("gray", "red"), ncol = 4)
dev.off()

pdf(paste0(out_dir, "known_marker/germ_featureplot.pdf"), width = 28, height = 7)
FeaturePlot(object = data.integrated, features = MARKERS$germ, min.cutoff = 0, max.cutoff = c(4, 3,  5), cols = c("gray", "red"), ncol = 4)
dev.off()

pdf(paste0(out_dir, "known_marker/endotherial_featureplot.pdf"), width = 28, height = 7)
FeaturePlot(object = data.integrated, features = MARKERS$endothelial, min.cutoff = 0, max.cutoff = c(4, 4, 4), cols = c("gray", "red"), ncol = 3)
dev.off()

pdf(paste0(out_dir, "known_marker/erythoro_featureplot.pdf"), width = 28, height = 7)
FeaturePlot(object = data.integrated, features = MARKERS$erythro, min.cutoff = 0, max.cutoff = c(10, 10, 3), cols = c("gray", "red"), ncol = 3)
dev.off()

pdf(paste0(out_dir, "known_marker/mk_featureplot.pdf"), width = 28, height = 7)
FeaturePlot(object = data.integrated, features = MARKERS$mk, min.cutoff = 0, max.cutoff = c(2, 2.5, 3.5), cols = c("gray", "red"), ncol = 3)
dev.off()

