
#==============================================================================================#

#intiial settings

#==============================================================================================#
# name pattern
run_id <- "E12_iD6"

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
    E12 = "./data/mapped/E12",
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
data_list[["E12"]] <- subset(data_list[["E12"]], subset = nFeature_RNA > 750 & nFeature_RNA < 6500 & nCount_RNA < 50000 & percent.mt < 25)
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
    replace(which(date == "E12"), 12) %>%
    replace(which(date == "iD6"), 12) -> age

vivo.vitro <- data.integrated@meta.data$sample
vivo.vitro %>%
    replace(which(vivo.vitro == "E12"), "vivo") %>%
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

# Clustering & UMAP (Supplementary Figure 6E)

#==============================================================================================#
#clustering
data.integrated <- RunPCA(object = data.integrated, npcs = 16, verbose = FALSE)
data.integrated <- RunUMAP(object = data.integrated, reduction = "pca", dims = 1:16)
data.integrated <- FindNeighbors(object = data.integrated, dims = 1:16)
data.integrated <- FindClusters(object = data.integrated, algorithm = 1, resolution = 0.2)

save(data.integrated, file = paste0(data_dir, "data.integrated.QCed.ccreg.clust.RData"))

# UMAP colored by cluster and splited by sample
pdf(paste0(out_dir, "by.sample_umap.pdf"), height=7, width=14)
DimPlot(object = data.integrated, reduction = "umap", split.by="sample", pt.size = 1.2)
dev.off()

#==============================================================================================#

#  Maeker gene expression (Supplementary Figure 6F)

#==============================================================================================#
dir.create(paste0(out_dir, "known_marker"))

#Germ cells
multiVlnPlot(object=data.integrated, features = MARKERS$germ, title="Germ Cell Markers", out_file =paste0(out_dir, "known_marker/germ_vln.pdf"))

#Endothelial
multiVlnPlot(object=data.integrated, features = MARKERS$endothelial, title="Endothelial Markers", out_file =paste0(out_dir, "known_marker/endothelial_Vln.pdf"))

#erythrovyte
multiVlnPlot(object=data.integrated, features =  MARKERS$erythro, title="Erythrocyte Markers",out_file =paste0(out_dir, "known_marker/ery_Vln.pdf"))

#megakaryocyte
multiVlnPlot(object=data.integrated, features = MARKERS$mk, title="MK Markers",out_file =paste0(out_dir, "known_marker/mk_vln.pdf"))

#==============================================================================================#

# Extraction of Gonadal Somatic Cells (Figure 2C, D)

#==============================================================================================#
data.integrated.sel <- subsetUmapClust(
    data =data.integrated,
    subset.name = "seurat_clusters",
    subset.value = 0:5,
    resolution=0.2,  # resolution shoud be optimized
    npcs = 20,
    dims = 1:20
)

save(data.integrated.sel, file = paste0(data_dir, "data.integrated.QCed.ccreg.clust.sel.RData"))

#UMAP colored by cluster and splited by sample
pdf(paste0(out_dir, "soma_by.sample_umap.pdf"), height=7, width=14)
DimPlot(object = data.integrated.sel, reduction = "umap", split.by="sample",pt.size = 1.2)
dev.off()

#UMAP colored by cell cycle phase and splited by sample
pdf(paste0(out_dir, "soma_cc.umap.pdf"), height=7, width=14)
DimPlot(object = data.integrated.sel, reduction = "umap", group.by="Phase", split.by="sample",pt.size = 1.2)
dev.off()

#==============================================================================================#

#Proportion of each cluster (Figure 2E)

#==============================================================================================#
#Reorder the clusters (1(EP), 5(EP), 0(STpro), 3(ST), 2(GR), 4(GR))
levels(data.integrated.sel@meta.data$seurat_cluster) <-  c(1, 5, 0, 3, 2, 4)
Idents(data.integrated.sel) <- data.integrated.sel@meta.data$seurat_cluster

#computtion of proportion
prop_cluster <- prop.table(x=table(data.integrated.sel@meta.data$seurat_cluster, data.integrated.sel@meta.data$sample), margin =2)

#color setting
ncluster <- length(levels(data.integrated.sel@meta.data$seurat_cluster))
clust_col = gg_color_hue(ncluster)
names(clust_col) = paste0("S", c(0:5))
clust_col <- clust_col[ c("S1", "S5", "S0", "S3", "S2", "S4")]

#plot
Cluster_name=factor(paste0("S", rownames(prop_cluster)), levels=c("S1", "S5", "S0", "S3", "S2", "S4"))    # To order, making a factor object for the cluster
df <- data.frame(Cluster=Cluster_name, D6=prop_cluster[,2], E12.5=prop_cluster[,1])
df2 <- tidyr::gather(df, key=Sample, value=Proportion, -Cluster, factor_key = TRUE)

##vidualization
theme <- theme(panel.background = element_blank(),    # initialization
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            plot.title = element_text(size = 20,hjust = 0.5),
            axis.text.x = element_text(size=30),
            axis.text.y = element_text(size=20),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks = element_line(size=0.5),
            axis.line = element_line(size=0.5),
            legend.text = element_text(size=20),
            legend.title =  element_text(size=20),
            plot.margin = unit(c(1,1,1,1),"line")
    )

g <- ggplot(df2, aes(x= Sample, y = Proportion, fill=Cluster))
        g <- g + geom_bar(stat = "identity")
        g <- g + scale_y_continuous(labels = percent)    # display as %
        g <- g + scale_fill_manual(values = clust_col, limits = c("S0", "S1", "S2","S3", "S4", "S5"))
        g <- g + theme

        g <- g + annotate("text",size = 10, x = 1, y = 0.07, label = "GR")
        g <- g + annotate("text",size = 10, x = 1, y = 0.19, label = "GR")
        g <- g + annotate("text",size = 10, x = 1, y = 0.34, label = "ST")
        g <- g + annotate("text",size = 10, x = 1, y = 0.56, label = "STpro")
        g <- g + annotate("text",size = 10, x = 1, y = 0.706, label = "EP")
        g <- g + annotate("text",size = 10, x = 1, y = 0.865, label = "EP")

        g <- g + annotate("text",size = 10, x = 2, y = 0.023, label = "GR")
        g <- g + annotate("text",size = 10, x = 2, y = 0.27, label = "GR")
        g <- g + annotate("text",size = 10, x = 2, y = 0.532, label = "ST")
        g <- g + annotate("text",size = 10, x = 2, y = 0.675, label = "STpro")
        g <- g + annotate("text",size = 10, x = 2, y = 0.815, label = "EP")
        g <- g + annotate("text",size = 10, x = 2, y = 0.935, label = "EP")
        
ggsave(filename = paste0(out_dir, "soma_proportion.pdf"), plot = g, width = 12, height = 14)

#==============================================================================================#

#Soma correlation matrix (Figure 2F)

#==============================================================================================#
#add a meta data
data.integrated.sel@meta.data$sample_clust <- paste0(data.integrated.sel@meta.data$sample, "-", data.integrated.sel@meta.data$seurat_clusters)
Idents(data.integrated.sel) <- data.integrated.sel@meta.data$sample_clust

#comutation of the average expression of a cluster of each sample
av.exp <- AverageExpression(data.integrated.sel)$integrated

#correlation coefficient
cor.exp <- cor(av.exp)

#hierarchical clustering
hc <- hclust(dist(cor.exp))

#color setting for R
breaks1 <- seq(0.96, 1, length=500)    #the range should be optimized
mycols1 <- colorRamp2(breaks1,colorRampPalette(brewer.pal(9, "GnBu"))(500)) # color for motif enrichment

#cluster color setting
data_set <- unique(data.integrated.sel@meta.data$sample_clust)

cluster <- sapply(strsplit(data_set, "-"), function(x) x[2])
names(cluster) <- data_set

ncluster=length(unique(cluster))
clust_col = gg_color_hue(ncluster)
names(clust_col) <- sort(unique(cluster))

#cluster annotation
ha_row = rowAnnotation(df = data.frame(cluster = cluster), col = list(cluster = clust_col), width = unit(1, "cm"), annotation_legend_param = list(title_gp = gpar(fontsize = 12), labels_gp = gpar(fontsize = 12)))
ha_col = HeatmapAnnotation(df = data.frame(cluster = cluster), col = list(cluster = clust_col), height = unit(1, "cm"), annotation_legend_param = list(title_gp = gpar(fontsize = 12), labels_gp = gpar(fontsize = 12)))

#correlation matrix
ht1 <- Heatmap(cor.exp,
                col = mycols1,
                name = "R", #title of legend, 
                row_names_gp = gpar(fontsize = 12),    # Text size for row names
                column_names_gp = gpar(fontsize = 12),    # Text size for column names
                column_names_rot = 90,
                heatmap_legend_param =  list(title_gp = gpar(fontsize = 12), labels_gp = gpar(fontsize = 12)),
                row_dend_width = unit(2, "cm"),   #row dendrogram width
                column_dend_height = unit(2, "cm"),   #column dendrogram height
                row_dend_side = "left",
                column_dend_side = "top",
                cluster_rows = hc,
                cluster_columns = hc,
                top_annotation = ha_col,
                left_annotation = ha_row,
                width = unit(14, "cm"),
                height = unit(14, "cm"),
                cell_fun = function(j, i, x, y, width, height, fill){
                        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "black", fill = NA))
                }
            )

draw(ht1, merge_legend = TRUE)

pdf(paste0(out_dir, "soma_correlationplot.pdf"), height = 9, width = 9)
draw(ht1, merge_legend = TRUE)
dev.off()

#==============================================================================================#

# Soma Marker expression (Supplementary Figure 6G)

#==============================================================================================#
granulosa_markers <- c("Kitl", "Inha", "Foxl2", "Runx1")
stromal_markers <- c("Wnt5a", "Pdgfra", "Tcf21", "Acta2")
progenitor_marker <- c("Sox11",  "Ecm1", "Nr5a1")

pdf(paste0(out_dir, "known_marker/soma_granulosa.pdf"), width=28, height=7)
FeaturePlot(object = data.integrated.sel, features = granulosa_markers , min.cutoff = 0, max.cutoff = c(2.0, 2.5, 1.5, 2), cols=c("gray", "red"), ncol=4)
dev.off()

pdf(paste0(out_dir, "known_marker/soma_stromal.pdf"), width=28, height=7)
FeaturePlot(object = data.integrated.sel, features = stromal_markers, min.cutoff = 0, max.cutoff = c(3.0, 2.0, 3, 3.0), cols=c("gray", "red"), ncol=4)
dev.off()

pdf(paste0(out_dir, "known_marker/soma_progenitor.pdf"), width=21, height=7)
FeaturePlot(object = data.integrated.sel, features = progenitor_marker, min.cutoff = c(0,0,0),  max.cutoff = c(3.5, 3.0, 1.0), cols=c("gray", "red"), ncol=4)
dev.off()

#==============================================================================================#

#Soma Cluster specific gene expression (Supplementary Figure 6H)

#==============================================================================================#
#identification of cluster specific genes. We used MAST algorism.
library(dplyr)
integrated.sel.markers <- FindAllMarkers(
object = data.integrated.sel,
only.pos = FALSE,
min.pct = 0.25,
logfc.threshold = 0.25,
test.use = "MAST"
)

save(integrated.sel.markers, file=paste0(data_dir, "soma_markers.RData"))

#Heatmap of top5 genes
integrated.sel.markers %>%
filter(avg_logFC >= 0) %>%
group_by(cluster) %>%
top_n(5, avg_logFC) -> top5.posi

pdf(paste0(out_dir, "soma_top5.posi.heatmap.pdf"), , height=7, width=7)
DoHeatmap(object = data.integrated.sel, features=top5.posi$gene, size=3, label=FALSE)
dev.off()

#==============================================================================================#

# differentially expression analysis  (Supplementary Figure 6I)

#==============================================================================================#
#parameters
control <- "E12"
treatment <- "id6"
outname <-  "E12"

seurat_clusters <- levels(unlist(data.integrated.sel[["seurat_clusters"]]))

#DEG analysis (MAST)
diff_genes.list <- list()
for (i in seq(length(seurat_clusters))){
    data_subset <- eval(parse(text=paste0("subset(data.integrated.sel, subset = seurat_clusters == ", seurat_clusters[i], ")")))    #Extraction of the target cluster
    Idents(object =data_subset) <- data_subset@meta.data$sample
    
    diff_genes <- try(
            Deg(seurat_object=data_subset, 
            control = control, 
            treatment = treatment, 
            outID = paste0(outname, "_C", seurat_clusters[i]),
            only.pos = FALSE,
            min.pct = 0.25,
            logfc.threshold = 0,
            test.use = "MAST",
            min.cells.feature = 1,
            min.cells.group = 1),
        silent = FALSE)    #DEG analysis
    if(class(diff_genes) == "try-error"){    # In case DEG analysis is error
        diff_genes <- NA
    }

    diff_genes.list[[as.numeric(i)]] <- diff_genes
    names(diff_genes.list)[as.numeric(i)] <- seurat_clusters[i]
}

#sig DEG list
sig_genes.list <- lapply(diff_genes.list, function(x){
    if(is.data.frame(x) && nrow(x) > 0){
        x %>%
            mutate(gene = rownames(x)) %>%
            dplyr::filter(avg_logFC >= 0.25 | avg_logFC <= -0.25) %>%
            dplyr::filter(p_val_adj < 10e-2)
    }
})

#MA plot
p1.list <- list()
for(i in seurat_clusters){
    sig_genes <- eval(parse(text=paste0("sig_genes.list$'", i, "'$gene")))
    if(length(sig_genes) == 0) sig_genes <- NULL

    title <- paste0("Cluster ", i)
    
    p1 <- averageDegMA(seurat_object=data.integrated.sel,
        cluster_key = "seurat_clusters", 
        clust = i, 
        comp_key = "sample", 
        control = control, 
        treatment =treatment,
        sig_genes = sig_genes,
        genes_to_label = NULL,
        title = title
    )
    listname <- paste0("Cluster_", i)
    p1.list <- eval(parse(text = paste0("c(p1.list, list(", listname,  " = p1))")))
}

#vidualization
p1.list2 <- list(p1.list$Cluster_0, p1.list$Cluster_1, p1.list$Cluster_5, p1.list$Cluster_3, p1.list$Cluster_2, p1.list$Cluster_4)
g1 <- plot_grid(plotlist = p1.list2, nrow = 1, align="h")
ggsave( paste0(out_dir, "soma_MA_all.pdf"), g1, height=7, width=42)