#we will re-analyse the scRNA-seq data GSE136831 
#download the data ---------
kaminski_new <- readRDS(file = "/home/isilon/users/o_syarif/COPD machine learning/Rdata/2024-02-16_scrna-v5.rds") #download the file 
##from GSE136831 https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE136831&format=file

#relabel the immune cells to match FACS cohort
kaminski_new$new_ID <- kaminski_new$Manuscript_Identity
kaminski_new$new_ID[which(kaminski_new$Manuscript_Identity=="ncMonocyte" )] <- "Non-classical monocytes"
kaminski_new$new_ID[which(kaminski_new$Manuscript_Identity=="cMonocyte" )] <- "Classical monocytes"
kaminski_new$new_ID[which(kaminski_new$Manuscript_Identity=="Macrophage_Alveolar" )] <- "Macrophage_Alveolar"
kaminski_new$new_ID[which(kaminski_new$Manuscript_Identity=="NK" )] <- "NK cells"
kaminski_new$new_ID[which(kaminski_new$Manuscript_Identity=="Macrophage" )] <- "Macrophage"
kaminski_new$new_ID[which(kaminski_new$Manuscript_Identity=="cDC2" )] <- "cDC2"
kaminski_new$new_ID[which(kaminski_new$Manuscript_Identity=="cDC1" )] <- "cDC1"
kaminski_new$new_ID[which(kaminski_new$Manuscript_Identity=="DC_Mature" )] <- "DC_Mature"
kaminski_new$new_ID[which(kaminski_new$Manuscript_Identity=="B" )] <- "B cells"
kaminski_new$new_ID[which(kaminski_new$Manuscript_Identity=="B_Plasma" )] <- "B Plasma"
kaminski_new$new_ID[which(kaminski_new$Manuscript_Identity=="T" )] <- "T helper"
kaminski_new$new_ID[which(kaminski_new$Manuscript_Identity=="T_Cytotoxic" )] <- "T cytotoxic"
kaminski_new$new_ID[which(kaminski_new$Manuscript_Identity=="T_Regulatory" )] <- "Tregs"
kaminski_new$new_ID[which(kaminski_new$Manuscript_Identity=="Mast" )] <- "Mast cells"
kaminski_new$new_ID[which(kaminski_new$Manuscript_Identity=="pDC" )] <- "pDC"
kaminski_new$new_ID[which(kaminski_new$Manuscript_Identity=="ILC_A" )] <- "ILC_A"
kaminski_new$new_ID[which(kaminski_new$Manuscript_Identity=="ILC_B" )] <- "ILC_B"
kaminski_new$new_ID[which(kaminski_new$Manuscript_Identity=="DC_Langerhans" )] <- "DC_Langerhans"


unique(kaminski_new$Disease_Identity)
#remove IPF samples
kaminski_new <- subset(kaminski_new, Disease_Identity=="IPF", invert=T)
#remove two samples with low cell count
kaminski_new <- subset(kaminski_new, Subject_Identity %in% c("137CO", "152CO"), invert=T)

#save file for downstream analysis
save(kaminski_new, file="/home/isilon/users/o_syarif/COPD machine learning/COPD_paper/Fig3H/input/scRNA_COPD_control.RData")

#preprocess the data ---------
#load the required libraries
library(RColorBrewer)
library(readxl)
library(tidyverse)

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

gc()
kaminski_new <- FindNeighbors(kaminski_new, reduction = "harmony", dims = 1:30)
kaminski_new <- RunUMAP(kaminski_new, reduction = "harmony", dims = 1:30, 
                        reduction.name = "umap.harmony")

library(ggplot2)
a<-DimPlot(
  kaminski_new,
  reduction = "umap.harmony",
  group.by = c("Manuscript_Identity"),label.size = 4, repel = T, 
  label = T, split.by="Disease_Identity")&NoAxes()+
  theme(legend.position = "none")+ggtitle(NULL);a

datainputdir <- "/home/isilon/users/o_syarif/COPD machine learning/COPD_paper/Fig3H/input/"
save(kaminski_new, file = paste0(datainputdir, "scRNA_COPD_control_preprocessed.RData"))

# Fig 3B UMAP ---------------------------------------------------
datainputdir <- "/home/isilon/users/o_syarif/COPD machine learning/COPD_paper/Fig3H/input/"
load(paste0(datainputdir, "scRNA_COPD_control_preprocessed.RData"))

#upload the new metadata for new label and colors
load("/home/isilon/users/o_syarif/COPD machine learning/COPD_paper/Fig5/revision_iScience/script/final revision/Fig3/input data/metadata.RData")
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
#colors for immune cells
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

#colors for structural cells
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
names(colors_dim) <- merge_color$num_clust

library(ggplot2)
a <- DimPlot(kaminski_new, group.by = "num_clust", split.by = "Disease_Identity", #change new_name to num_clust to get the number label
             label = T,cols = colors_dim, 
             repel = T)&NoAxes()&ggtitle("")&guides(color=guide_legend(ncol=4));a
a_plot<- a+theme(legend.position = "none")
library(ggplotify)
library(ggpubr)
legend <- get_legend(a+theme(legend.position = "right"))
legend <- as_ggplot(legend)
legend
ggsave(file=paste0(as.character(Sys.Date()),"_", "umap_scrna.png"), plot=a_plot,
       device = "png", units = c("in"), width = 7, height = 4, dpi = 300)
ggsave(file=paste0(as.character(Sys.Date()),"_", "umap_scrna LG.png"), plot=legend,
       device = "png", units = c("in"), width = 10, height = 4, dpi = 300)

#Fig 3C PCA biplot -------------------------
library(RColorBrewer)
library(readxl)
library(tidyverse)

#load the frequency data
datainputdir <- "/home/isilon/users/o_syarif/COPD machine learning/COPD_paper/Fig2/input/"
#upload the immune cell frequency table 
scRNA <- paste0(datainputdir, "2024-02-16_cell_freq.xlsx")
cell.freq <- read_excel(scRNA, sheet = 1)

annotation_col <- cell.freq %>% 
  filter(Disease_Identity != "IPF") %>% # remove IPF samples
  filter(!Subject_Identity %in% c("137CO", "152CO")) %>% #remove the two samples with low cell counts
  select(Disease_Identity) 
annotation_col <- as.data.frame(annotation_col)
annotation_col$Disease_Identity<-as.factor(annotation_col$Disease_Identity)

rownames(annotation_col) <- paste0("a_",  c(1:44))

scrna.df <- cell.freq  %>%  
  filter(Disease_Identity != "IPF") %>%
  filter(!Subject_Identity %in% c("137CO", "152CO")) %>%
  select(-Disease_Identity, -Subject_Identity)

scrna.df <- log(scrna.df+1)

rownames(scrna.df) <- paste0("a_",  c(1:44))

#annotation
color <- c( "darkgrey", "darkred")
names(color) <-  c("Control","COPD")


#callback = function(hc, ...){dendsort(hc)}
#rename the annotation
colnames(scrna.df) <- case_when(colnames(scrna.df) =="B" ~ "B cells",
                                colnames(scrna.df) =="B_Plasma"  ~ "B Plasma",
                                colnames(scrna.df) =="ILC_A" ~ "ILC_A",
                                colnames(scrna.df) =="ILC_B" ~ "ILC_B",
                                colnames(scrna.df) =="Macrophage" ~ "Macrophage",
                                colnames(scrna.df) =="Macrophage_Alveolar" ~ "Macrophage_Alveolar",
                                colnames(scrna.df) =="NK" ~ "NK cells",
                                colnames(scrna.df) =="T" ~ "T helper",
                                colnames(scrna.df) =="T_Cytotoxic" ~ "T cytotoxic",
                                colnames(scrna.df) =="T_Regulatory" ~ "Tregs",
                                colnames(scrna.df) =="cDC2" ~ "cDC2",
                                colnames(scrna.df) =="cMonocyte" ~ "Classical monocytes",
                                colnames(scrna.df) =="ncMonocyte" ~ "Non-classical monocytes",
                                colnames(scrna.df) =="pDC" ~ "pDC",
                                colnames(scrna.df) =="DC_Langerhans" ~ "DC_Langerhans",
                                colnames(scrna.df) =="Mast" ~ "Mast cells",
                                colnames(scrna.df) =="cDC1" ~ "cDC1",
                                colnames(scrna.df) =="DC_Mature" ~ "DC_Mature")

#calculate PCA
PCA_modelobject <- prcomp(scrna.df, scale=TRUE, center=TRUE)
summary(PCA_modelobject)


theme_journal <- theme(axis.line.x = element_line(color="black", size = 0.5),
                       axis.line.y = element_line(color="black", size = 0.5))+
  theme(axis.text.x = element_text(colour = "black", size = 14, angle = 0, hjust=1, vjust=1))+
  theme(axis.text.y = element_text(colour = "black", size = 14))+
  theme(axis.title = element_text(size = 14, colour = "black"))+
  theme(legend.position="right")+
  theme(legend.text = element_text(size = 14))+
  theme(legend.title = element_text(size = 16))

plotdata_cPCA$Disease_Identity <- factor(plotdata_cPCA$Disease_Identity, levels = c("Control","COPD" ))
pov <- PCA_modelobject$sdev^2/sum(PCA_modelobject$sdev^2)
barplot(pov)
#visualisation
library(ggplot2)
theme_set(theme_pubr(base_size=6, border = T, legend = "none"))
theme_journal <- theme(panel.grid.major = element_line(color = "light grey", linewidth = 0.2, linetype = "solid"),
                       panel.grid.minor = element_line(color = "light grey", linewidth = 0.05, linetype = "solid"),
                       plot.title = element_text(size = rel(1)),
                       axis.line.x = element_line(color="black", linewidth = 0.2),
                       axis.line.y = element_line(color="black", linewidth = 0.2),
                       axis.text = element_text(color="black", size = 6),
                       panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.2)
)
library("factoextra")


#note biplot works only for cPCA and mPCA but not mixO onces (some kind of class error in the delivered PCA objects)
Biplot_cPCA <- fviz_pca_biplot(PCA_modelobject, title = "PCA-Biplot", 
                               habillage= plotdata_cPCA$Disease_Identity,
                               mean.point = FALSE, pointsize = 1, pointshape = 19,
                               label ="var", select.var = list(contrib = 10), #addEllipses=TRUE, ellipse.level=0.95
                               col.var = "black",labelsize = 2, arrowsize = 0.2, repel = TRUE) +  #can't color by contrib since uses already color_scale for PCA scores
  scale_color_manual(values=color) +
  expand_limits(x = c(min(PCA_modelobject$x[,"PC1"])*1.1, max(PCA_modelobject$x[,"PC1"])*1.1),
                y = c(min(PCA_modelobject$x[,"PC2"])*1.1, max(PCA_modelobject$x[,"PC2"])*1.1)) + 
  guides(color = guide_legend(override.aes = list(size = 1))) + #reduces dot size in legend
  xlab(paste0("PC1 ", round(pov[1]*100, 1), "%"))+ylab(paste0("PC2 ", round(pov[2]*100,1), "%"))+ 
  theme_journal +   # +coord_fixed() #makes quadratic
  theme(legend.position="none", plot.title = element_text(color="black", size = 6),
        plot.caption = element_text(size = 4, hjust = 0), plot.caption.position = "plot",
        aspect.ratio = 1,
        axis.title.x = element_text(colour = "black", size = 6),
        axis.title.y = element_text(colour = "black", size = 6) );Biplot_cPCA


ggsave(paste0("scrna__cPCA_biplot.png"),
       plot = Biplot_cPCA , 
       device = png(), width=65, height=70, units = "mm", dpi = 900)
dev.off()

#Fig 3D Violin Plot ------------------

#Fig 3E Dot Plot --------------------
datainputdir <- "/home/isilon/users/o_syarif/COPD machine learning/COPD_paper/Fig3H/input/"
load(paste0(datainputdir, "scRNA_COPD_control_preprocessed.RData"))
library(Seurat)
kaminski_new <- JoinLayers(kaminski_new)
detected_cytokines <- c("CCL2", "CCL3", "CCL4","CCL5", "CCL11", "CCL17", "CCL19", "CCL20", "CXCL1", "CXCL5", "CXCL9", "CXCL10",
                        "CXCL11", "CXCL12", "CXCL13", "CSF2", "IFNL1", "IFNB1", "IFNG","IL1B", "IL6", "CXCL8","IL10",
                        "TNF", "LTA","TSLP", "LTB")

#we calculate the significance of each cytokine using loop function
markers <- data.frame(matrix(NA,    # Create empty data frame
                             nrow = 1,
                             ncol = 7))
colnames(markers) <- c("p_val","avg_log2FC", "pct.1","pct.2","p_val_adj","population", "gene")

for (i in unique(kaminski_new$new_ID)){
  tryCatch({
    Idents(kaminski_new) <- kaminski_new$new_ID
    teta <- subset(kaminski_new, idents = c(i))
    teta <- NormalizeData(teta)
    Idents(teta) <- teta$Disease_Identity
    obj.markers <- FindMarkers(teta, ident.1 = c("COPD"), ident.2 = c("Control"), 
                               features = detected_cytokines, method="wilcox")
    obj.markers$population <- i
    obj.markers$gene <- rownames(obj.markers)
    markers <- rbind(markers, obj.markers)
  }, error = function(e) {
    # Code to handle the error (skip the iteration)
    message(paste("Error occurred at i =", i, "; skipping..."))
  })
}

markers <- markers[-1,]
markers_wilcox <- markers

#No TSLP pass min.pct threshold

x <- markers_wilcox$p_val_adj[which(markers_wilcox$p_val_adj >0)]
dat <- markers_wilcox %>%
  mutate(p_val_adj = ifelse(p_val_adj == 0, min(x), p_val_adj)) # we change the 0 p value to be minimum values of p value 
dat$gene <- factor(dat$gene, levels = detected_cytokines)

#now we move to the FACS dataset
#FACS data
datainputdir <- "/home/isilon/users/o_syarif/COPD machine learning/COPD_paper/Fig3H/input/"
data_FACS <- read.csv(file = paste0(datainputdir, "FACS_foldchange.csv"))
data_FACS$X <- NULL

# add additional line for the dotplot wilcoxon graph
#after matching with data_FACS data we ended up having the following matched cytokines between plasma, lung, and scRNA
not_overlap <- c("CCL19_L", "CXCL13_L", "IFN_la2_3_S" ,"IFNa2_S",  "IFNb_S", "TNF_C_L")
overlap_cytokine <- c("CCL11_L","CCL11_S",
                      "CCL17_L","CCL17_S",                   
                      "CCL2_L","CCL2_S",
                      "CCL20_L","CCL20_S",                   
                      "CCL3_L","CCL3_S",
                      "CCL4_L","CCL4_S",                    
                      "CCL5_L","CCL5_S",                     
                      "CXCL1_L","CXCL1_S",
                      "CXCL10_L","CXCL10_S",
                      "CXCL11_L","CXCL11_S",
                      "CXCL12_L","CXCL12_S",                   
                      "CXCL5_L","CXCL5_S",
                      "CXCL9_L","CXCL9_S", 
                      "GM_CSF_L","GM_CSF_S",
                      "IFN_la1_L","IFN_la1_S",
                      "IFNg_L","IFNg_S",                    
                      "IL10_L","IL10_S",
                      "IL1b_L","IL1b_S",
                      "IL6_L","IL6_S",
                      "IL8_L","IL8_S",
                      "TNFa_L","TNFa_S",
                      "TNFb_L","TNFb_S")

data <- data_FACS[which(data_FACS$.y. %in% overlap_cytokine),]
temp <- as.data.frame(str_split_fixed(data$.y., "_", 2))
data$cytokine <- str_remove_all(data$.y., "_L|_S")
cytokines_FACS_available <- overlap_cytokine

data <- data[which(data$.y. %in% cytokines_FACS_available),]
data$origin <- ifelse(grepl("_S", data$.y.), "plasma", "lung")
data$cytokine <- str_replace_all(data$cytokine,"IFNg", "IFNG")
data$cytokine <- str_replace_all(data$cytokine,"IL1b", "IL1B")
data$cytokine <- str_replace_all(data$cytokine,"TNFa", "TNF")
data$cytokine <- str_replace_all(data$cytokine,"TNFb", "LTA")
data$cytokine <- str_replace_all(data$cytokine,"IFN_la1", "IFNL1")
data$cytokine <- str_replace_all(data$cytokine,"GM_CSF", "CSF2")
data$cytokine <- str_replace_all(data$cytokine,"IL8", "CXCL8")

#to check again the availability of the cyokines in scRNA seq
unique(dat$gene)

intercept <- unique(intersect(unique(dat$gene), unique(data$cytokine)))

cytokine_plot <- intercept
data_toplot <- data[, c("p","p.adj","origin","log2_FC", "cytokine")]
data_toplot$origin <- ifelse(data_toplot$origin=="plasma", "Plasma", "Lung")

#take the scRNA data
dat_toplot <- dat[, c("p_val","p_val_adj","population","avg_log2FC", "gene")]
colnames(dat_toplot)
colnames(dat_toplot) <- c("p","p.adj","origin","log2_FC", "cytokine")
sort(unique(dat_toplot$origin))

#continue FACS data
data_toplot <- data_toplot[which(data_toplot$cytokine %in% cytokine_plot),]
dat_toplot <- dat_toplot[which(dat_toplot$cytokine %in% cytokine_plot),]

#combine FACS and scRNA data
data_bind <- rbind(dat_toplot, data_toplot)

#filter based on the p values
x <- data_bind$p.adj[which(data_bind$p.adj >0)]
datg <- data_bind %>%
  mutate(p.adj = ifelse(p.adj == 0, min(x), p.adj)) # we change the 0 p value to be minimum values of p value 

datg$padjust <- -log10(datg$p.adj)
datg$padjust[datg$padjust >20] <- 20

unique(kaminski_new$new_ID)
order_immune <- c("B cells", "B Plasma", "cDC1", "cDC2",
                  "DC_Langerhans", "DC_Mature",
                  "ILC_A", "ILC_B", "Macrophage", 
                  "Macrophage_Alveolar", "Mast cells","NK cells", "Non-classical monocytes", "Classical monocytes",
                  "pDC", "T cytotoxic", "T helper",
                  "Tregs")
order_structural <- c("ATI","ATII",
                      "Basal", "Ciliated","Club",
                      "Fibroblast","Goblet",
                      "Lymphatic" , 
                      "Mesothelial",         "Myofibroblast",
                      "Pericyte","PNEC","SMC",
                      "VE_Arterial","VE_Capillary_A"  ,
                      "VE_Capillary_B","VE_Peribronchial","VE_Venous")

datg$origin <- factor(datg$origin, levels = c(order_structural, order_immune, "Lung", "Plasma"))
datg$cytokine <- factor(datg$cytokine, levels = c("CCL2" ,  "CCL3",   "CCL4" ,  "CCL5"  , "CCL11" ,
                                                  "CCL17",  "CCL20" , "CXCL1", "CXCL5", "CXCL8", 
                                                  "CXCL9" , "CXCL10", "CXCL11",
                                                  "CXCL12", "CSF2",   "IFNL1"  ,"IFNG" ,  "IL1B",   "IL6",
                                                  "IL10"  , "TNF"   , "LTA"))

#plot both data
library(RColorBrewer)
color_palette <- c(rev(colorRampPalette(brewer.pal(9, "Blues"))(48)), "#ffffff", "#ffffff" , 
                   colorRampPalette(brewer.pal(9, "Reds"))(48))
min_boundary <- min(datg$log2_FC,na.rm=TRUE)
max_boundary <- max(datg$log2_FC,na.rm=TRUE)
color_breaks <- c(seq(min_boundary, min_boundary/49, -min_boundary/49), 
                  0, 
                  seq(max_boundary/49, max_boundary, max_boundary/49)) #asymmetric sclae maxing out contrast for up and down
color_palette <- c(rev(colorRampPalette(brewer.pal(9, "Blues"))(48)), "#ffffff", "#ffffff" , 
                   colorRampPalette(brewer.pal(9, "Reds"))(48))


#datg <- na.omit(datg)
dd <- datg[which(datg$origin %in% c("Plasma", "Lung")),]
limit <- max(abs(dd$log2_FC)) * c(-1, 1)

#plot
dd$source <- "cytokine cohort"
c <- ggplot(dd, aes(y=cytokine, x=origin)) +
  geom_point(aes(size=padjust, color =log2_FC))+ 
  scale_color_gradientn(colours = color_palette, limits=limit) +
  theme_classic()+coord_flip()+theme_bw(base_size = 16)+
  facet_wrap(vars(source),ncol=1, scales = "free_y")+
  theme(strip.text.x = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 11))+ 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "gray",linetype = "dashed", size=0.35), 
        panel.border = element_rect(colour = "black", fill=NA, size=2))+
  scale_size_continuous(range = c(-1,6), breaks = c(0, 10, 20), labels = c("0", "10", ">20"))

c_plot <- c+theme(legend.position = "none")
c_plot
library(ggplotify)
library(ggpubr)
a <- get_legend(c)
a <- as_ggplot(a)
a
ggsave(filename = paste0(plotdir, Sys.Date(),"_","dotplots_GSE136831_COPDcontrol_cytokine_wilcox_lung_plasma",".png"),
       device = "png", width = 6, plot = c_plot, dpi=300,
       height =1.3,
       units = c("in"))
ggsave(filename = paste0(plotdir, Sys.Date(),"_","dotplots_GSE136831_COPDcontrol_cytokine_wilcox_lung_plasma_legend",".png"),
       device = "png", width = 3, plot = a, dpi=300,
       height =5,
       units = c("in"))

dd <- datg[-which(datg$origin %in% c("Plasma", "Lung")),]
limit <- max(abs(dd$log2_FC)) * c(-1, 1)
library(forcats)

dd$source <- ifelse(dd$origin %in% order_immune, "immune", "structural")
c <- ggplot(dd, aes(y=cytokine, x=fct_rev(origin))) +
  geom_point(aes(size=padjust, color =log2_FC))+ 
  scale_color_gradientn(colours = color_palette, limits=limit) +
  theme_classic()+coord_flip()+theme_bw(base_size = 16)+
  facet_wrap(vars(source),ncol=1, scales = "free_y")+
  theme(strip.text.x = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))+ 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "gray",linetype = "dashed", size=0.35), 
        panel.border = element_rect(colour = "black", fill=NA, size=2))+
  scale_size_continuous(range = c(-1,6), breaks = c(0, 10, 20), labels = c("0", "10", ">20"));c
c_plot <- c+theme(legend.position = "none")
c_plot
a <- get_legend(c)
a <- as_ggplot(a)
a
ggsave(filename = paste0(plotdir, Sys.Date(),"_","dotplots_GSE136831_COPDcontrol_cytokine_wilcox_onlyscRNA",".png"),
       device = "png", width = 6.5, plot = c_plot, dpi=300,
       height =9,
       units = c("in"))

ggsave(filename = paste0(plotdir, Sys.Date(),"_","dotplots_GSE136831_COPDcontrol_cytokine_wilcox_onlyscRNA_legend",".png"),
       device = "png", width = 3, plot = a, dpi=300,
       height =5,
       units = c("in"))



#Fig 3F interaction plot -------------
#receptor ligand analysis
library(dplyr)
library(readxl)
library(Seurat)

#cytokine of interest
cytokine <- c("CCL5","CXCL9","CXCL10")

#upload data required for the cell frquency
datainputdir <- "/home/isilon/users/o_syarif/COPD machine learning/COPD_paper/Fig2/input/"

scRNA <- paste0(datainputdir, "2024-02-16_cell_freq.xlsx")
cell.freq <- read_excel(scRNA, sheet = 1)

#annotation
color <- c( "darkgrey", "darkred")
names(color) <-  c("Control","COPD")



all_cytokines <- c("CCL5","CXCL9","CXCL10")

library(tidyverse)
a <- cell.freq %>% pivot_longer(where(is.numeric)) %>%
  filter(!Subject_Identity %in% c("137CO", "152CO")) %>%
  filter(!Disease_Identity %in% c("IPF")) %>%
  mutate(Disease_Identity = factor(Disease_Identity, levels = c("COPD", "Control"))) %>%
  mutate(value = log(value+1))
a
dg <-a %>% group_by(name, Disease_Identity) %>%
  summarise(across(where(is.numeric), median))

#rename immune cells annotation to match FACS cohort
dg$name <- case_when(dg$name =="B" ~ "B cells",
                     dg$name =="B_Plasma"  ~ "B Plasma",
                     dg$name =="ILC_A" ~ "ILC_A",
                     dg$name =="ILC_B" ~ "ILC_B",
                     dg$name =="Macrophage" ~ "Macrophage",
                     dg$name =="Macrophage_Alveolar" ~ "Macrophage_Alveolar",
                     dg$name =="NK" ~ "NK cells",
                     dg$name =="T" ~ "T helper",
                     dg$name =="T_Cytotoxic" ~ "T cytotoxic",
                     dg$name =="T_Regulatory" ~ "Tregs",
                     dg$name =="cDC2" ~ "cDC2",
                     dg$name =="cMonocyte" ~ "Classical monocytes",
                     dg$name =="ncMonocyte" ~ "Non-classical monocytes",
                     dg$name =="pDC" ~ "pDC",
                     dg$name =="DC_Langerhans" ~ "DC_Langerhans",
                     dg$name =="Mast" ~ "Mast cells",
                     dg$name =="cDC1" ~ "cDC1",
                     dg$name=="DC_Mature" ~ "DC_Mature")


#upload ramilowski LR database 
ramilwoski <- read.csv(file = "/home/isilon/users/o_syarif/COPD machine learning/COPD_paper/Fig5/revision_iScience/ramilowski_LR.csv")
ramilowski_sub <- ramilwoski[which(ramilwoski$ligand %in% all_cytokines),]
ramilowski_sub <-unique(ramilowski_sub)

#cell cell interaction calculation
library(Seurat)
datainputdir <- "/home/isilon/users/o_syarif/COPD machine learning/COPD_paper/Fig3H/input/"
load(paste0(datainputdir, "scRNA_COPD_control.RData"))
Idents(kaminski_new) <- kaminski_new$new_ID
sbset_cell <- c("T cytotoxic", "T helper", "Macrophage", 
                "Macrophage_Alveolar", "Classical monocytes", 
                "Non-classical monocytes")
kaminski_new <- subset(kaminski_new, new_ID %in% sbset_cell)
selection_ctrl <- subset(kaminski_new, Disease_Identity == "Control") #change to COPD or Control
selection_copd <- subset(kaminski_new, Disease_Identity == "COPD") #change to COPD or Control

selection <- selection_ctrl
DefaultAssay(selection) <- "RNA"
selection <- NormalizeData(selection)
df_all <- as.matrix(GetAssayData(selection, slot="data", assay = "RNA"))
genes <- as.character(unique(c(ramilowski_sub$ligand, ramilowski_sub$receptor)))
df_all <- df_all[which(rownames(df_all) %in% genes),]
df_all <- as.data.frame(t(df_all))
df_all <- cbind(df_all, selection$new_ID)
colnames(df_all)[13] <- "cellannottaion"
#used minimum 1% of expressing cell type to define the interactions
perc=1 
perc=perc/100
result=data.frame()
res=data.frame()
i=1
j=2

df_all$cellannottaion <- as.factor(df_all$cellannottaion)

#loop over each cluster to find pairs
for(i in 1:(length(levels(df_all$cellannottaion)))){
  for(j in 1:(length(levels(df_all$cellannottaion )))){
    #from the large martix, subselect receptor and lig subgoups (if i=1 and j=2, keep cells in grps 1 and 2)
    test=df_all[df_all$cellannottaion==levels(df_all$cellannottaion)[i] | df_all$cellannottaion==levels(df_all$cellannottaion)[j],]
    #Subselect genes in receptor list in cells in rec subgroup (say 1)
    R_c1=test[test$cellannottaion==levels(df_all$cellannottaion)[i] ,(colnames(test) %in% ramilowski_sub$receptor)]
    #Subselect genes in ligand list in cells in lig subgroup (say 2)
    L_c2=test[test$cellannottaion==levels(df_all$cellannottaion)[j] , (colnames(test) %in% ramilowski_sub$ligand)]
    if(nrow(R_c1)!=0 &nrow(L_c2)!=0){
      #keep genes that are expressed in more than user-input percent of the cells
      keep1 = colnames(R_c1)[which(colSums(R_c1>0)>=perc*dim(R_c1)[1])]
      keep2 = colnames(L_c2)[which(colSums(L_c2>0)>=perc*dim(L_c2)[1])]
      R_c1=R_c1[,keep1]
      L_c2=L_c2[,keep2]
      #get list of lig-rec pairs
      res=ramilowski_sub[(ramilowski_sub$ligand %in% colnames(L_c2)) & (ramilowski_sub$receptor %in% colnames(R_c1)),]
    }else{}
    if(nrow(res)!=0){
      res$Receptor_cluster=levels(df_all$cellannottaion)[i]
      res$Lig_cluster=levels(df_all$cellannottaion)[j]
      result=rbind(result,res)
    }else{result=result}
  }
}

results_COPD <- result
unique(dg$Disease_Identity)
dg_COPD <- dg[which(dg$Disease_Identity=="Control"),]
dg_COPD <- dg_COPD[which(dg_COPD$name %in% sbset_cell),]
dg_COPD$Receptor_cluster <- dg_COPD$name
merge_COPD <- left_join(results_COPD, dg_COPD, by="Receptor_cluster")
merge_COPD$flow <-1

library(ggalluvial)
merge_COPD$pair <-paste0(merge_COPD$ligand,"_", merge_COPD$receptor)
dy <- merge_COPD[, c("Disease_Identity","value","Receptor_cluster","Lig_cluster","pair")]
dy$celltype <-paste0(dy$Receptor_cluster,"/", dy$Lig_cluster)
dy$Receptor_cluster <-NULL
dy$Lig_cluster <-NULL

dy <- separate_rows(dy, celltype, sep="/")

#pseudobulk of receptor and ligand expression
pseudo <-AggregateExpression(selection_ctrl, assays = "RNA", group.by = "new_ID", return.seurat = T)
matgene <- GetAssayData(pseudo, assay="RNA",slot="data")
matgene <- as.matrix(matgene)
matgene <- matgene[unique(merge_COPD$ligand),]
library(reshape2)
matgene <- melt(matgene)
matgene$Var2 <- as.character(matgene$Var2)
matgene$Var2[which(matgene$Var2=="Macrophage-Alveolar")] <- "Macrophage_Alveolar"

merge_COPD$gg <- paste0(merge_COPD$Lig_cluster, "/", merge_COPD$ligand)
matgene$gg <- paste0(matgene$Var2, "/", matgene$Var1)
matgene$Var1 <-NULL
matgene$Var2 <-NULL
colnames(matgene)[1] <-"weight"

pseudo <-AggregateExpression(selection_ctrl, assays = "RNA", group.by = "new_ID", return.seurat = T)
matgener <- GetAssayData(pseudo, assay="RNA",slot="data")
matgener <- as.matrix(matgener)
matgener <- matgener[unique(merge_COPD$receptor),]
matgener <- melt(matgener)
matgener$Var2 <- as.character(matgener$Var2)
matgener$Var2[which(matgener$Var2=="Macrophage-Alveolar")] <- "Macrophage_Alveolar"


merge_COPD$ggx <- paste0(merge_COPD$Receptor_cluster, "/", merge_COPD$receptor)
matgener$ggx <- paste0(matgener$Var2, "/", matgener$Var1)
matgener$Var1 <-NULL
matgener$Var2 <-NULL


colnames(matgener)[1] <-"weight_receptor"
merge_xr <- left_join(merge_COPD, matgener, by="ggx")
merge_xg <- left_join(merge_xr, matgene, by="gg")
colnames(merge_xg)

#prepare for table of ligand and receptor
dt_copd <- merge_xg[, c("ligand"       ,    "receptor"    ,     "Receptor_cluster" ,"Lig_cluster" , #change dt_ctrl for control data
                        "Disease_Identity", "value" , "weight_receptor",  "weight" )]
colnames(dt_copd) <- c("ligand"       ,    "receptor"    ,     "receptor_cluster" ,"ligand_cluster" ,
                       "diagnosis", "receptor_cluster_proportion" , "receptor_expression",  "ligand_expression" )

dk <- dg_COPD
colnames(dk)[4] <- "ligand_cluster"
dt_copd <- left_join(dt_copd, dk[, c(3,4)], by="ligand_cluster")
colnames(dt_copd)[9] <- "ligand_cluster_proportion"

#save the data
library(openxlsx)
wa <- createWorkbook()
addWorksheet(wa, "RL_ctrl")
addWorksheet(wa, "RL_copd")


# Write data to the worksheet (example: write a data frame)
writeData(wa, sheet = "RL_ctrl", dt_ctrl)
writeData(wa, sheet = "RL_copd", dt_copd)

saveWorkbook(wa, file = paste0("/home/isilon/users/o_syarif/COPD machine learning/COPD_paper/Fig5/revision_iScience/output/RLdata.xlsx"), overwrite = TRUE)

#end saving 

merge_xg <- merge_xg[, c("Lig_cluster", "Receptor_cluster","weight", "weight_receptor")]

merge_xg$score <- merge_xg$weight*merge_xg$weight_receptor

#visualisation
library(igraph)
merge_xg_sub <- merge_xg[, c("Lig_cluster", "Receptor_cluster", "score")]
colnames(merge_xg_sub) <- c("from","to","weight")
df_copd <- merge_xg_sub
merge_xg_sub$weight <- merge_xg_sub$weight*40


#copd
set.seed(123)
g <- graph_from_data_frame(merge_xg_sub, directed = T)
cells <-names(V(g))
dg_COPD$name <- factor(dg_COPD$name , levels = cells)
d <- dg_COPD[order(dg_COPD$name),]

scores <- d$value

V(g)$weight <- scores

V(g)$color <- c("#96DED1","#ADD8E6","#A7C7E7","#9FE2BF","#FA5F55","#FAA0A0")
macrophage ="#ADD8E6"
mac_alv ="#A7C7E7"
cmono ="#96DED1"
ncmono= "#9FE2BF"
Thelper = 	"#FAA0A0"
Tcyto ="#FA5F55"

#we only choose certain populations as the top deregulated cell types based on calculations on Fig 3E
populations_to_focus <- c("Macrophage"     ,         "Macrophage_Alveolar" ,
            "Non-classical monocytes", "Classical monocytes",
            "T cytotoxic"    ,        "T helper")
# Plot the weighted graph
png("graph_COPD2.png", width = 10, height = 10, res=300, units = c("in"))
par(las = 2, mar = c(2, 2, 10, 1)) 
b <-plot(g,
         edge.width = E(g)$weight, # Edge width based on weight
         vertex.size =  V(g)$weight*20,
         edge.arrow.size=0.2,
         vertex.label=NA,
         vertex.label.color = "black",
         edge.color = "grey30",
         layout = layout_in_circle(g,order = populations_to_focus)
)
dev.off()


#ctrl 
set.seed(123)
g <- graph_from_data_frame(merge_xg_sub, directed = T)
cells <-names(V(g))
dg_COPD$name <- factor(dg_COPD$name , levels = cells)
d <- dg_COPD[order(dg_COPD$name),]

scores <- d$value

V(g)$weight <- scores
V(g) #please check the color code for each celltype
V(g)$color <- c("#ADD8E6","#A7C7E7","#9FE2BF","#96DED1","#FA5F55","#FAA0A0")
#macrophage ="#ADD8E6"
#mac_alv ="#A7C7E7"
#cmono ="#96DED1"
#ncmono= "#9FE2BF"
#T helper = 	"#FAA0A0"
#Tcyto ="#FA5F55"

populations_to_focus <- c("Macrophage"     ,         "Macrophage_Alveolar" ,
            "Non-classical monocytes", "Classical monocytes",
            "T cytotoxic"    ,        "T helper")

# Plot the weighted graph
png("graph_control2.png", width = 10, height = 10, res=300, units = c("in"))
par(las = 2, mar = c(2, 2, 10, 1)) 
b <-plot(g,
         edge.width = E(g)$weight, # Edge width based on weight
         vertex.size =  V(g)$weight*20,
         edge.arrow.size=0.2,
         #vertex.label=NA,
         vertex.label.color = "black",
         edge.color = "grey30",
         layout = layout_in_circle(g,order = populations_to_focus)
)
dev.off()

