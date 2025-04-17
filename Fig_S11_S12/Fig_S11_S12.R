#scRNA--------
#load the data------------------
library(ggplot2)
library(dplyr)
library(stringr)
library(Seurat)

#take the non filtered data containing only COPD and controls see script Fig3_scRNA.R
load("non_integrated_preprocessessed_scrna_GSE136831.RData")

##non integration------------
kaminski_new <- NormalizeData(kaminski_new)
kaminski_new <- FindVariableFeatures(kaminski_new)
kaminski_new <- ScaleData(kaminski_new)
kaminski_new <- RunPCA(kaminski_new)
kaminski_new <- FindNeighbors(kaminski_new, reduction = "pca", dims = 1:30)
kaminski_new <- RunUMAP(kaminski_new, reduction = "pca", dims = 1:30, 
                        reduction.name = "umap.pca")



###metadata ---------------
dir <- paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/")
load(file =paste0(dir, "input data/metadata.RData"))
names(colorsss)[38] <-"Ionocyte"
df_order$new_name <- paste0("(", df_order$order_num, ")", df_order$species)

library(Seurat)
order_immune <- c("B cells", "B Plasma", "cDC1", "cDC2",
                  "Classical monocytes", "DC_Langerhans", "DC_Mature",
                  "ILC_A", "ILC_B", "Macrophage", 
                  "Macrophage_Alveolar", "Mast cells","NK cells", "Non-classical monocytes",
                  "pDC", "T cytotoxic", "T helper",
                  "Tregs")
order_structural <- c("ATI","ATII","Aberrant_Basaloid" ,"Ionocyte" ,
                      "Basal", "Ciliated","Club",
                      "Fibroblast","Goblet",
                      "Lymphatic" , 
                      "Mesothelial",         "Myofibroblast",
                      "Pericyte","PNEC","SMC",
                      "VE_Arterial","VE_Capillary_A"  ,
                      "VE_Capillary_B","VE_Peribronchial","VE_Venous")

color_df<- as.data.frame(colorsss)
color_df$species <- rownames(color_df)
merge_color <- merge(color_df, df_order, by="species")

colors_dim <- merge_color$colorsss
names(colors_dim) <- merge_color$new_name

merge_color$species <- factor(merge_color$species, levels = c(order_immune, order_structural))
merge_color <- merge_color[order(merge_color$species),]
level_dim <-unique(merge_color$new_name)
merge_color$num_clust <- as.character(merge_color$order_num)
kaminski_new$new_name <- merge_color$new_name[match(kaminski_new$new_ID,
                                                    merge_color$species)]
kaminski_new$num_clust <- merge_color$num_clust[match(kaminski_new$new_ID,
                                                      merge_color$species)]
kaminski_new$new_name <- factor(kaminski_new$new_name, levels = level_dim)

merge_color$new_name
merge_color$new_color <- NA
merge_color$new_color[which(merge_color$new_name=="(1)B cells")] <- "#E3735E"
merge_color$new_color[which(merge_color$new_name=="(2)B Plasma")] <- "#A42A04"
merge_color$new_color[which(merge_color$new_name=="(3)cDC1" )] <- "#C1E1C1"
merge_color$new_color[which(merge_color$new_name=="(4)cDC2")] <- "#93C572"
merge_color$new_color[which(merge_color$new_name=="(5)Classical monocytes")] <-"#96DED1"
merge_color$new_color[which(merge_color$new_name=="(6)DC_Langerhans" )] <- "#8A9A5B"
merge_color$new_color[which(merge_color$new_name=="(7)DC_Mature" )] <- "#40B5AD"
merge_color$new_color[which(merge_color$new_name=="(8)ILC_A" )] <- "#CCCCFF"
merge_color$new_color[which(merge_color$new_name=="(9)ILC_B")] <- "#C3B1E1"
merge_color$new_color[which(merge_color$new_name=="(10)Macrophage" )] <- "#ADD8E6"
merge_color$new_color[which(merge_color$new_name=="(11)Macrophage_Alveolar")] <- "#A7C7E7"
merge_color$new_color[which(merge_color$new_name=="(12)Mast cells")] <- "#BDB5D5"
merge_color$new_color[which(merge_color$new_name=="(13)NK cells")] <- "#6082B6"
merge_color$new_color[which(merge_color$new_name=="(14)Non-classical monocytes")] <- "#9FE2BF"
merge_color$new_color[which(merge_color$new_name=="(15)pDC" )] <- "#708090"
merge_color$new_color[which(merge_color$new_name=="(16)T cytotoxic")] <- "#FA5F55"
merge_color$new_color[which(merge_color$new_name=="(17)T helper")] <- "#FAA0A0"
merge_color$new_color[which(merge_color$new_name== "(18)Tregs")] <- "#F88379"

merge_color$new_color[which(merge_color$new_name=="(19)ATI")] <-  "#F8C8DC"
merge_color$new_color[which(merge_color$new_name=="(20)ATII")] <-"#F8C8DC"
merge_color$new_color[which(merge_color$new_name=="(21)Aberrant_Basaloid")] <-  "#F8C8DC"
merge_color$new_color[which(merge_color$new_name=="(22)Ionocyte")] <-"#F8C8DC"
merge_color$new_color[which(merge_color$new_name=="(23)Basal" )] <- "#F8C8DC"
merge_color$new_color[which(merge_color$new_name=="(24)Ciliated" )] <- "#F8C8DC"
merge_color$new_color[which(merge_color$new_name=="(25)Club")] <- "#F8C8DC"
merge_color$new_color[which(merge_color$new_name=="(26)Fibroblast")] <- "#FAD5A5"
merge_color$new_color[which(merge_color$new_name=="(27)Goblet")] <- "#F8C8DC"
merge_color$new_color[which(merge_color$new_name=="(28)Lymphatic")] <- "#C2B280"
merge_color$new_color[which(merge_color$new_name=="(29)Mesothelial")] <- "#FAD5A5"
merge_color$new_color[which(merge_color$new_name=="(30)Myofibroblast")] <- "#FAD5A5"
merge_color$new_color[which(merge_color$new_name=="(31)Pericyte")] <- "#FAD5A5"
merge_color$new_color[which(merge_color$new_name=="(32)PNEC")] <- "#F8C8DC"
merge_color$new_color[which(merge_color$new_name=="(33)SMC" )] <- "#FAD5A5"
merge_color$new_color[which(merge_color$new_name=="(34)VE_Arterial")] <- "#C2B280"
merge_color$new_color[which(merge_color$new_name=="(35)VE_Capillary_A")] <- "#C2B280"
merge_color$new_color[which(merge_color$new_name=="(36)VE_Capillary_B" )] <- "#C2B280"
merge_color$new_color[which(merge_color$new_name=="(37)VE_Peribronchial")] <- "#C2B280"
merge_color$new_color[which(merge_color$new_name=="(38)VE_Venous")] <- "#C2B280"

### UMAP non integrated---------------------------------------------------
colors_dim <- merge_color$new_color
names(colors_dim) <- merge_color$new_name
#names(colors_dim) <- merge_color$num_clust

library(ggplot2)
a<-DimPlot(
  kaminski_new,
  reduction = "umap.pca",
  group.by = c("Subject_Identity"),label.size = 4, repel = F, 
  label = F)&NoAxes()+
  theme(legend.position = "right")+ggtitle(NULL)
a_plot <- a+theme(legend.position = "none")
ggsave(filename = paste0(Sys.Date(),"_","UMAp_scrna_subject_pre",".png"),
       device = "png", width = 5, plot = a_plot, dpi=300,
       height =5,
       units = c("in"))
ggsave(filename = paste0(Sys.Date(),"_","UMAp_scrna_subject_pre LG",".png"),
       device = "png", width = 9, plot = a, dpi=300,
       height =5,
       units = c("in"))


a<-DimPlot(
  kaminski_new,
  reduction = "umap.pca",cols = colors_dim, 
  group.by = c("new_name"),label.size = 4, repel = F, 
  label = F)&NoAxes()+
  theme(legend.position = "right")+ggtitle(NULL)
a_plot <- a+theme(legend.position = "none")
ggsave(filename = paste0(Sys.Date(),"_","UMAp_scrna_celltype_pre",".png"),
       device = "png", width = 5, plot = a_plot, dpi=300,
       height =5,
       units = c("in"))
ggsave(filename = paste0(Sys.Date(),"_","UMAp_scrna_celltype_pre LG",".png"),
       device = "png", width = 9, plot = a, dpi=300,
       height =5,
       units = c("in"))

a<-DimPlot(
  kaminski_new,
  reduction = "umap.pca",
  group.by = c("Disease_Identity"),label.size = 4, repel = F,
  cols = c( "darkgrey", "darkred"), 
  label = F)&NoAxes()+
  theme(legend.position = "right")+ggtitle(NULL)
a_plot <- a+theme(legend.position = "none")
ggsave(filename = paste0(Sys.Date(),"_","UMAp_scrna_diagnosis_pre",".png"),
       device = "png", width = 5, plot = a_plot, dpi=300,
       height =5,
       units = c("in"))
ggsave(filename = paste0(Sys.Date(),"_","UMAp_scrna_diagnosis_pre LG",".png"),
       device = "png", width = 9, plot = a, dpi=300,
       height =5,
       units = c("in"))

plotoutput <- "/home/isilon/users/o_syarif/COPD machine learning/COPD_paper/data_output/"

save(kaminski_new, file = paste0(plotoutput, "non_integrated_preprocessessed_scrna_GSE136831.RData"))

##integration -----------
#take the filtered data containing only COPD and controls see script Fig3_scRNA.R
load("scRNA_COPD_control.RData")

kaminski_new[["RNA"]] <- split(kaminski_new[["RNA"]], f = kaminski_new$Subject_Identity)
library(Seurat)
kaminski_new <- NormalizeData(kaminski_new)
kaminski_new <- FindVariableFeatures(kaminski_new)
kaminski_new <- ScaleData(kaminski_new)
kaminski_new <- RunPCA(kaminski_new)
kaminski_new <- IntegrateLayers(
  object = kaminski_new, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = T
)

### UMAP integrated---------------------------------------------------
dir <- paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/")
load(file =paste0(dir, "input data/metadata.RData"))
names(colorsss)[38] <-"Ionocyte"
df_order$new_name <- paste0("(", df_order$order_num, ")", df_order$species)

library(Seurat)
order_immune <- c("B cells", "B Plasma", "cDC1", "cDC2",
                  "Classical monocytes", "DC_Langerhans", "DC_Mature",
                  "ILC_A", "ILC_B", "Macrophage", 
                  "Macrophage_Alveolar", "Mast cells","NK cells", "Non-classical monocytes",
                  "pDC", "T cytotoxic", "T helper",
                  "Tregs")
order_structural <- c("ATI","ATII","Aberrant_Basaloid" ,"Ionocyte" ,
                      "Basal", "Ciliated","Club",
                      "Fibroblast","Goblet",
                      "Lymphatic" , 
                      "Mesothelial",         "Myofibroblast",
                      "Pericyte","PNEC","SMC",
                      "VE_Arterial","VE_Capillary_A"  ,
                      "VE_Capillary_B","VE_Peribronchial","VE_Venous")
#kaminski_new$new_label <- ifelse(kaminski_new$new_ID %in% order_immune, kaminski_new$new_ID, "structural")
#kaminski_new$new_label <- factor(kaminski_new$new_label, levels = c(order_immune, "structural"))

color_df<- as.data.frame(colorsss)
color_df$species <- rownames(color_df)
merge_color <- merge(color_df, df_order, by="species")

colors_dim <- merge_color$colorsss
names(colors_dim) <- merge_color$new_name

merge_color$species <- factor(merge_color$species, levels = c(order_immune, order_structural))
merge_color <- merge_color[order(merge_color$species),]
level_dim <-unique(merge_color$new_name)
merge_color$num_clust <- as.character(merge_color$order_num)
kaminski_new$new_name <- merge_color$new_name[match(kaminski_new$new_ID,
                                                    merge_color$species)]
kaminski_new$num_clust <- merge_color$num_clust[match(kaminski_new$new_ID,
                                                      merge_color$species)]
kaminski_new$new_name <- factor(kaminski_new$new_name, levels = level_dim)

merge_color$new_name
merge_color$new_color <- NA
merge_color$new_color[which(merge_color$new_name=="(1)B cells")] <- "#E3735E"
merge_color$new_color[which(merge_color$new_name=="(2)B Plasma")] <- "#A42A04"
merge_color$new_color[which(merge_color$new_name=="(3)cDC1" )] <- "#C1E1C1"
merge_color$new_color[which(merge_color$new_name=="(4)cDC2")] <- "#93C572"
merge_color$new_color[which(merge_color$new_name=="(5)Classical monocytes")] <-"#96DED1"
merge_color$new_color[which(merge_color$new_name=="(6)DC_Langerhans" )] <- "#8A9A5B"
merge_color$new_color[which(merge_color$new_name=="(7)DC_Mature" )] <- "#40B5AD"
merge_color$new_color[which(merge_color$new_name=="(8)ILC_A" )] <- "#CCCCFF"
merge_color$new_color[which(merge_color$new_name=="(9)ILC_B")] <- "#C3B1E1"
merge_color$new_color[which(merge_color$new_name=="(10)Macrophage" )] <- "#ADD8E6"
merge_color$new_color[which(merge_color$new_name=="(11)Macrophage_Alveolar")] <- "#A7C7E7"
merge_color$new_color[which(merge_color$new_name=="(12)Mast cells")] <- "#BDB5D5"
merge_color$new_color[which(merge_color$new_name=="(13)NK cells")] <- "#6082B6"
merge_color$new_color[which(merge_color$new_name=="(14)Non-classical monocytes")] <- "#9FE2BF"
merge_color$new_color[which(merge_color$new_name=="(15)pDC" )] <- "#708090"
merge_color$new_color[which(merge_color$new_name=="(16)T cytotoxic")] <- "#FA5F55"
merge_color$new_color[which(merge_color$new_name=="(17)T helper")] <- "#FAA0A0"
merge_color$new_color[which(merge_color$new_name== "(18)Tregs")] <- "#F88379"

merge_color$new_color[which(merge_color$new_name=="(19)ATI")] <-  "#F8C8DC"
merge_color$new_color[which(merge_color$new_name=="(20)ATII")] <-"#F8C8DC"
merge_color$new_color[which(merge_color$new_name=="(21)Aberrant_Basaloid")] <-  "#F8C8DC"
merge_color$new_color[which(merge_color$new_name=="(22)Ionocyte")] <-"#F8C8DC"
merge_color$new_color[which(merge_color$new_name=="(23)Basal" )] <- "#F8C8DC"
merge_color$new_color[which(merge_color$new_name=="(24)Ciliated" )] <- "#F8C8DC"
merge_color$new_color[which(merge_color$new_name=="(25)Club")] <- "#F8C8DC"
merge_color$new_color[which(merge_color$new_name=="(26)Fibroblast")] <- "#FAD5A5"
merge_color$new_color[which(merge_color$new_name=="(27)Goblet")] <- "#F8C8DC"
merge_color$new_color[which(merge_color$new_name=="(28)Lymphatic")] <- "#C2B280"
merge_color$new_color[which(merge_color$new_name=="(29)Mesothelial")] <- "#FAD5A5"
merge_color$new_color[which(merge_color$new_name=="(30)Myofibroblast")] <- "#FAD5A5"
merge_color$new_color[which(merge_color$new_name=="(31)Pericyte")] <- "#FAD5A5"
merge_color$new_color[which(merge_color$new_name=="(32)PNEC")] <- "#F8C8DC"
merge_color$new_color[which(merge_color$new_name=="(33)SMC" )] <- "#FAD5A5"
merge_color$new_color[which(merge_color$new_name=="(34)VE_Arterial")] <- "#C2B280"
merge_color$new_color[which(merge_color$new_name=="(35)VE_Capillary_A")] <- "#C2B280"
merge_color$new_color[which(merge_color$new_name=="(36)VE_Capillary_B" )] <- "#C2B280"
merge_color$new_color[which(merge_color$new_name=="(37)VE_Peribronchial")] <- "#C2B280"
merge_color$new_color[which(merge_color$new_name=="(38)VE_Venous")] <- "#C2B280"


colors_dim <- merge_color$new_color
names(colors_dim) <- merge_color$new_name


library(ggplot2)
a<-DimPlot(
  kaminski_new,
  group.by = c("Subject_Identity"),label.size = 4, repel = F, 
  label = F)&NoAxes()+
  theme(legend.position = "right")+ggtitle(NULL)
a_plot <- a+theme(legend.position = "none")
ggsave(filename = paste0(Sys.Date(),"_","UMAp_scrna_subject_post",".png"),
       device = "png", width = 5, plot = a_plot, dpi=300,
       height =5,
       units = c("in"))
ggsave(filename = paste0(Sys.Date(),"_","UMAp_scrna_subject_post LG",".png"),
       device = "png", width = 9, plot = a, dpi=300,
       height =5,
       units = c("in"))


a<-DimPlot(
  kaminski_new,
  cols = colors_dim, 
  group.by = c("new_name"),label.size = 4, repel = F, 
  label = F)&NoAxes()+
  theme(legend.position = "right")+ggtitle(NULL)
a_plot <- a+theme(legend.position = "none")
ggsave(filename = paste0(Sys.Date(),"_","UMAp_scrna_celltype_post",".png"),
       device = "png", width = 5, plot = a_plot, dpi=300,
       height =5,
       units = c("in"))
ggsave(filename = paste0(Sys.Date(),"_","UMAp_scrna_celltype_post LG",".png"),
       device = "png", width = 9, plot = a, dpi=300,
       height =5,
       units = c("in"))

a<-DimPlot(
  kaminski_new,
  group.by = c("Disease_Identity"),label.size = 4, repel = F,
  cols = c( "darkgrey", "darkred"), 
  label = F)&NoAxes()+
  theme(legend.position = "right")+ggtitle(NULL)
a_plot <- a+theme(legend.position = "none")
ggsave(filename = paste0(Sys.Date(),"_","UMAp_scrna_diagnosis_post",".png"),
       device = "png", width = 5, plot = a_plot, dpi=300,
       height =5,
       units = c("in"))
ggsave(filename = paste0(Sys.Date(),"_","UMAp_scrna_diagnosis_post LG",".png"),
       device = "png", width = 9, plot = a, dpi=300,
       height =5,
       units = c("in"))

#Spatial transcriptomics -------------
library(GeoMXAnalysisWorkflow)
library(standR)
library(GeomxTools)
library(Seurat)
library(SpatialDecon)
library(patchwork)
library(SpatialExperiment)
library(SingleCellExperiment)

#1) load the data
dir <- paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/")
load(file =paste0(dir, "input data/spatial_experiment_parenchyma.RData"))

#2) then continue with log transformation for pre batch correction
spe@assays@data@listData[["log_batch"]] <- log2(spe@assays@data@listData[["batch_corrected"]]) #change count to batch_corrected to check pre and post batch correction

#3) calculate and plot PCA
library(SpatialExperiment)
set.seed(100)
spe <- scater::runPCA(spe, exprs_values = "log_batch")
pca_results_tmm <- reducedDim(spe, "PCA")
spe@colData@listData[["batch"]] <- spe@colData@listData[["batch"]]
spe@colData@listData[["batch_fixation"]]  <- paste0(spe@colData@listData[["fixation"]], "/", spe@colData@listData[["batch"]] )
plotPairPCA(spe, precomputed = pca_results_tmm, color = location) 
a <-plotPairPCA(spe, precomputed = pca_results_tmm, color = batch)
a
a <-plotPairPCA(spe, precomputed = pca_results_tmm, color = fixation)
a
a <-plotPairPCA(spe, precomputed = pca_results_tmm, color = gold_broad)
a

#sve PCA
library(ggplotify)
plot <-as.ggplot(a);plot
ggsave(file=paste0(as.character(Sys.Date()),"_", "ST_batchcorrected_gold.png"), 
       plot=plot,
       device = "png", units = c("in"), width = 6, height = 6)



