#--------------------------------------------
# scCountProc
#--------------------------------------------
#' Reading cellranger processed 10x data as a seurat object.
#' 
#' Read the cell ranger processed 10x data, assign a sample name, normalize, and find veriable genes
#' 
#' @param data_path path of ceranger preocessed result (one upper from "outs" directory)
#' @param sample_name charector a sample name
#' @param normalization logical, if TRUE, normalization is performed
#' @param nfetures integer, number of variable genes to be found
#' 
#' @importFrom Seurat Read10X CreateSeuratObject RenameCells NormalizeData FindVariableFeatures
#' 
#' @return a seurat object
#' 
#' @keywords seurat 10x
#' @export 

scCountProc <- function(data_path, sample_name = sample_name, normalization = FALSE, nfeatures =2000){
    ## Count 
    library(Seurat)
    data.counts <- Read10X(data.dir = paste0(data_path, "/outs/filtered_feature_bc_matrix/"))
    data <- CreateSeuratObject(counts = data.counts)    # conversion to Seurat object
    data$sample <- sample_name    #サンプル名の付加
    data <- RenameCells(object=data, add.cell.id = sample_name)
    if (normalization == TRUE){
        data <- NormalizeData(object = data)    #正規化
    }
    data <- FindVariableFeatures(object = data, selection.method = "vst", nfeatures = nfeatures)    #変動の大きい遺伝子の抽出
    return(data)
}

#--------------------------------------------
# QC
#--------------------------------------------
#' QC for 10x data
#' 
#' this function output ggplot object to vidualize nFeature_RNA, nCount_RNA, and percent of mitochondria genes
#' 
#' @param data a Seurat objects
#' 
#' @importFrom Seurat PercentageFeatureSet VlnPlot FeatureScatter CombinePlots
#' @importFrom ggplot2 ggsave
#' 
#' @return a list of data of percent.mt and ggplot objects
#' 
#' @keywords seurat 10x
#' @export 

QC <- function(data){
    data[["percent.mt"]] <- PercentageFeatureSet(object = data, pattern = "^mt-")
    vp <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    sp <- CombinePlots(plots = list(
        FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt"),
        FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    ))
    out <- list(percent.mt = data[["percent.mt"]], violin_plot = vp, scatter_plot = sp)
    return(out)
}

#--------------------------------------------
# ccReg
#--------------------------------------------
#' Cell cycle regression
#' 
#' this function regress the effect of cell cycle
#' 
#' @param data a Seurat object
#' 
#' @importFrom Seurat PercentageFeatureSet VlnPlot FeatureScatter CombinePlots
#' @importFrom ggplot2 ggsave
#' 
#' @return a list of data of percent.mt and ggplot objects
#' 
#' @keywords seurat 10x
#' @export 

ccReg <- function(data){
    #load cell cycle markers(Tirosh et al, 2015)
    s.genes <- sapply(cc.genes.updated.2019$s.genes, function(x){ paste0(substr(x, 1, 1), tolower(substr(x, 2, nchar(x)))) })
    g2m.genes <- sapply(cc.genes.updated.2019$g2m.genes, function(x){ paste0(substr(x, 1, 1), tolower(substr(x, 2, nchar(x)))) })
    
    #Initialization step
    cat("Initialization step running...\n")
    data <- ScaleData(object = data, verbose = FALSE)
    data <- RunPCA(object = data, npcs = 30, verbose = FALSE)
    
    #cell-cycle scoring
    cat("cell-cycle scoring running...\n")
    data <- CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    
    #PCA before regression
    cat("PCA before regression...\n")
    data <- RunPCA(data, features = c(s.genes, g2m.genes))
    
    pca1 <- DimPlot(data)
    ggsave(pca1, file="ccReg_results/CC.pca.pdf")
    
    #Cell-Cycle regression
    cat("Cell Cycle Regression running...\n")
    data$CC.Difference <- data$S.Score - data$G2M.Score
    data <- ScaleData(data, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(data))
    # data <- ScaleData(data, vars.to.regress = "CC.Difference", features = rownames(data))
    
    #PCA after regression
    cat("PCA after regression...\n")
    data <- RunPCA(data, features = VariableFeatures(data), nfeatures.print = 10)
    data <- RunPCA(data, features = c(s.genes, g2m.genes))
    
    pca2 <- DimPlot(data)
    ggsave(pca2, file="ccReg_results/reg_CC.pca.pdf")
    
    cat("cell-cycle regression succeeded! \n")
    cat("all output files are stored in ccReg_results \n")
    return(data)
}

#--------------------------------------------
# subsetUmapClust
#--------------------------------------------
#' Extract subset and do clustering
#' 
#' Extract the subsets from a seurat object and do scaling, and PCA, UMAP and tSNE dimensionality reductions
#' 
#' @param data A seurat object
#' @param subset.name Parameter to subset on. Eg, the name of a gene, PC_1, a column name in object@meta.data, etc. Any argument that can be retreived using FetchData
#' @param subset.value Returns cells with the subset name equal to this value
#' @param resolution resolution for clustering
#' 
#' @importFrom Seurat Read10X CreateSeuratObject RenameCells NormalizeData FindVariableFeatures
#' 
#' @return a seurat object
#' 
#' @keywords seurat 10x
#' @export 
#' 
subsetUmapClust <- function(data,  subset.name = "seurat_clusters", subset.value, npcs = 30, dims = 1:30, resolution = 0.6){
    if(length(subset.value) == 1){
        data.subset <- eval(parse(text=paste0("subset(data, subset = ", subset.name, " == ", subset.value, ")" )))
    }else{
        multisub <- paste0("subset(data, ", 
                                paste(paste(subset.name, paste0("'",subset.value,"'"), sep="=="), collapse = "|"),
                                ")"
        )
        data.subset <- eval(parse(text=multisub))
    }
    data.subset <- RunPCA(object = data.subset, npcs = npcs, verbose = FALSE)
    data.subset <- RunUMAP(object = data.subset, reduction = "pca", dims = dims)
    data.subset <- FindNeighbors(object = data.subset, dims = dims)
    data.subset <- FindClusters(object = data.subset, reduction.type = "umap", resolution =  resolution)
    return(data.subset)
}

#--------------------------------------------
# gg_color_hue
#--------------------------------------------
#' Get the Default ggplot2 Color Palette
#' 
#' \code{gg_color_hue} returns the default ggplot color palette for a given
#' number of colors
#'
#' @author John Colby \url{http://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette}
#
#' @param n Number of Coolors
#' 
#' @export

gg_color_hue <- function(n){
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

#--------------------------------------------
# Deg 
#--------------------------------------------
#' compute defferentially expressed genes and output as a text file
#' 
#' compute defferentially expressed genes and output as a text file
#' 
#' @param seurat_object a seurat object
#' @param control control for the comarison. a element of the comp key.
#' @param treatment comparison for the comarison. a element of the comp key.
#' @param outID ID for output files

#' 
#' @importFrom Seurat FindMarkers 
#' 
#' @return a dataframe of differentially expressed genes
#' 
#' @keywords DEG
#' @export 

Deg <- function(seurat_object, control, treatment, outID,...){
    diff_genes <- FindMarkers(object=seurat_object, ident.1 = treatment, ident.2 = control,...)

    out_file <- paste0(outID, "_diff_genes_", treatment, ".", control, ".txt")
    write.table(diff_genes, file=out_file, sep="\t", quote=FALSE)

    return(diff_genes)
}

#--------------------------------------------
# averageDegMA
#--------------------------------------------
#' plot a MA-plot for two average gene expression data sets
#' 
#'  plot a MA-plot for two average gene expression data sets.
#' 
#' @param seurat_object a seurat object
#' @param cluster_key key for categorization. a header of meta.data of the seurat object
#' @param clust a charactor of cluster
#' @param comp_key key for comparison. a header of meta.data of the seurat object
#' @param control control for the comarison. a element of the comp key.
#' @param treatment comparison for the comarison. a element of the comp key.
#' @param sig_genes = a vector of significant gene names. Points correspondng to these genes are painted as red.
#' @param genes_to_label a vector of gene symbols to be labeled.
#' @param title a charactor. title of the MA-plot.

#' 
#' @importFrom Seurat Idents RotatedAxis DotPlot AverageExpression LabelPoints
#' @importFrom ggplot2 element_blank theme ggtitle xlim ylim geom_point ggplot
#' @importFrom cowplot plot_grid
#' @importFrom dplyr mutate %>%
#' 
#' @return a pdf file of MA-plot
#' 
#' @keywords DEG, MA-plot
#' @export 

averageDegMA <- function(seurat_object, 
    cluster_key = "seurat_clusters", 
    clust,
    comp_key, 
    control, 
    treatment, 
    sig_genes = NULL, 
    genes_to_label = NULL,
    title = NULL
    ){
    library(ggplot2)
    library(cowplot)
    message("extracting a subet data...")
    data_subset <- eval(parse(text=paste0("subset(seurat_object, subset = ", cluster_key, " == clust)")))
    Idents(data_subset) <- eval(parse(text=paste0("data_subset@meta.data$", comp_key)))

    #at least one sample does not have any cells, skip the comparison
    sample_included <- eval(parse(text=paste0(
        "all(c(any(data_subset@meta.data$",
        comp_key,
        " == treatment",
        "), any(data_subset@meta.data$",
        comp_key,
        " == control",
        ")))"
    )))

    if(sample_included == FALSE){
        warning("The cluster does not have the cell.", immediate. = TRUE)
        p1 <- NULL
    }else{
        #create averaged data model
        message("creating the averaged data model...")
        avg.data_subset <- log1p(AverageExpression(data_subset, verbose = FALSE)$RNA)
        avg.data_subset <- as.data.frame(apply(avg.data_subset, 2, function(x){replace(x, which(x==0), 1e-06)}))

        R <- round(cor(avg.data_subset[,treatment],avg.data_subset[,control]), digits = 2)
        nsig <- length(sig_genes)

        M <- eval(parse(text=paste0("avg.data_subset$", treatment, "- avg.data_subset$", control)))
        A <- apply(avg.data_subset[,c(treatment, control)], 1, mean)
        MA <- data.frame(M,A)

        sigs <- factor(rep("Insignificant", nrow(MA)), levels = c("Insignificant", "Significant"))
        if(!is.null(sig_genes)){
            sigs <- replace(sigs, which(rownames(MA) %in% sig_genes), "Significant")
        }
        MA <- cbind(MA, sigs)
        MA <- MA[sort.list(MA$sigs),]

        #plot
        message("plotting the MA-Plot...")
        p1 <- ggplot(MA, aes(x = A, y = M, color=sigs))
        p1 <- p1 + geom_point(size=1.80, alpha = 0.4)
        p1 <- p1 + scale_color_manual(values = c(Insignificant="gray40", Significant="red"))
        if(is.null(title)){
            p1 <- p1 + ggtitle(paste0("Cluster ", clust, "_", control, " vs. ", treatment))
        }else{
            p1 <- p1 + ggtitle(title)
        }
        p1 <- p1 + geom_hline(yintercept=0, linetype="dashed", colour="gray30", size=1)
        p1 <- p1 + theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border = element_rect(colour="gray1", fill=NA),
                axis.text = element_text(size=18),
				axis.title = element_text(size=18),
				plot.title = element_text(size=18, face="bold", hjust = 0.5),
				aspect.ratio = 1,
				legend.text = element_text(size =18),
				legend.title = element_blank(),
                legend.key = element_blank(),
                legend.position = c(.18,.92),
                legend.background = element_rect(fill = "gray97")
                )
        if(length(genes_to_label) !=0){
            p1 <- LabelPoints(plot = p1, points = genes_to_label, repel = TRUE, colour="black", size=5)
        }
        p1 <- p1 + annotate("text", size = 6, x=Inf,y=-Inf, hjust=1.26,vjust=-3.5, label = paste0("R = ", R))
        p1 <- p1 + annotate("text", size = 6, x=Inf,y=-Inf, hjust=1.1,vjust=-1.5, label = paste0("Number of DEGs: ",as.character(nsig)) , colour="red")

        plot(p1)
        
        plot_file=paste0("Cluster ", clust, "_", control, ".", treatment, "_MAPlot.pdf")
        ggsave(file=plot_file, plot=p1, width=7, height=7)

        message("Done\n")
    }
    return(p1)
}

#--------------------------------------------
# multiVlnPlot
#--------------------------------------------
#' plot a integrated multiple violin plots from seurat object
#' 
#' plot a integrated multiple violin plots from seurat object. wrapper of VlnPlot
#' 
#' @param object Seurat object
#' @param features Features to plot (gene expression, metrics, PC scores, anything that can be retreived by FetchData)
#' @param pt.size Point size for geom_violin
#' @param log plot the feature axis on log scale
#' @param title a charactor. title of the plot
#' @param out_file a charactor. file name of exprted plot

#' 
#' @importFrom Seurat VlnPlot 
#' @importFrom ggplot2 scale_y_continuous scale_fill_discrete ylab theme element_blank element_text element_rect get_legend guides guide_legend ggdraw draw_label
#' @importFrom cowplot plot_grid
#' 
#' @return a pdf file of violon plots
#' 
#' @keywords VlnPlot, Violin plot
#' @export 

multiVlnPlot <- function(object, features, pt.size = 0, log = FALSE, title = NULL, out_file = NULL, ...){
    library(Seurat)
    library(ggplot2)
    library(cowplot)

    scaleFUN <- function(x) sprintf("%.1f", x)

    gplot.list <- list()
    for(i in 1:length(features)){
        vg <- VlnPlot(object = object, features = features[i], pt.size = pt.size, log = log,...)
        vg <- vg + scale_y_continuous(labels=scaleFUN)
        vg <- vg + scale_fill_discrete(name="Cluster")
        vg <- vg + ylab(features[i])
        vg <- vg + theme(legend.position = 'none',
                        axis.title.x = element_blank(),
                        axis.title.y = element_text(face="bold", vjust = 0.5),
                        plot.title = element_blank(),
                        axis.text.x = element_blank(),
                        panel.border = element_rect(colour = "gray1", fill = NA),
                        )
        gplot.list <- c(gplot.list, list(vg))
    }

    legend <- get_legend(
        gplot.list[[1]] +
            theme(legend.position = "bottom", legend.justification="center") +
            guides(fill = guide_legend(nrow = 1, label.position = "bottom"))
    )

    # now add the title
    if(!is.null(title)){
        gtitle <- ggdraw() + 
            draw_label(
                title,
                fontface = 'bold',
                x = 0.5,
                hjust = 0.5
            ) 
    }else gtitle <- NULL

    gplot_all <- plot_grid(plotlist = gplot.list, ncol=1, align = "v", axis = "tblr")
    gplot_all <- plot_grid(
        gtitle, 
        gplot_all,
        legend,
        nrow = 3, 
        rel_heights = c(0.5 / (length(features) + 0.5), length(features) / (length(features) + 0.5), 0.5 / (length(features) + 0.5))
        )

    plot(gplot_all)

    if (!is.null(out_file)){
        ggsave(file=out_file, plot=gplot_all, width=14, height=length(features) * 1.3 +0.5)
    }
    return (gplot_all)
}

