#==============================================================================================#

# load functions

#==============================================================================================#
source("R/functions.r")

#==============================================================================================#

# load libraries

#==============================================================================================#
#scRNAseq
library(Seurat)

#general
library(tidyverse)
library(dplyr)

#ggplot
library(ggplot2)
library(ggsci)
library(cowplot)

#heatmap
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

#==============================================================================================#

# narker genes 

#==============================================================================================#
MARKERS <- list(
    pan = "Nr5a1",
    granulosa = c("Foxl2", "Inha", "Kitl"),
    progenitor = c("Sox11", "Ecm1", "Nr2f1"),
    stromal = c("Wnt5a", "Pdgfra", "Tcf21"),
    germ = c("Pou5f1", "Ddx4", "Dppa5a"),
    endothelial = c("Pecam1", "Flt1", "Ecscr"),
    erythro = c("Hbb-y","Hba-a1", "Bpgm"),
    mk = c("Ppbp", "Plek", "Pf4")
)
#==============================================================================================#

# ggplot theme

#==============================================================================================#
theme <- theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    plot.margin = unit(c(1, 1, 1, 1), "line")
)
