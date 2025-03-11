#we want to create the figure 2 in paper
library(Seurat)
plotdir <- "/home/isilon/users/o_syarif/COPD machine learning/COPD_paper/Fig2/output/plots/"
datainputdir <- "/home/isilon/users/o_syarif/COPD machine learning/COPD_paper/Fig2/input/"
dataoutputdir <- "/home/isilon/users/o_syarif/COPD machine learning/COPD_paper/Fig2/output/data/"


#for publication we start from here!
#for visualisation we load the relabelled scRNA
load(paste0(datainputdir, "scrna_relabel.RData"))

## pheatmap ---------
#load the required libraries
library(pheatmap)
library(RColorBrewer)
library(readxl)
library(tidyverse)

#load the frequency data
scRNA <- paste0(datainputdir, "2024-02-16_cell_freq.xlsx")
cell.freq <- read_excel(scRNA, sheet = 1)

annotation_col <- cell.freq %>% 
  filter(Disease_Identity != "IPF") %>%
  filter(!Subject_Identity %in% c("137CO", "152CO")) %>%
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

color_palette <- c(rev(colorRampPalette(brewer.pal(9, "Blues"))(48)), "#ffffff", "#ffffff" , 
                   colorRampPalette(brewer.pal(9, "Reds"))(48))	
#annotation
color <- c( "darkgrey", "darkred")
names(color) <-  c("Control","COPD")
annotation_colors <- list(color)
names(annotation_colors) <- "Disease_Identity"

callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}


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
a <-pheatmap(t(scrna.df), cluster_rows = T, cluster_cols = T, scale = "row", color = color_palette, 
             annotation_colors = annotation_colors,
             clustering_callback = callback,
             annotation_col = annotation_col,
             cellheight = 30, cellwidth = 10, fontsize = 12,
             show_colnames = F, show_rownames = T,
             #labels_row =  parse.t(temp.t$expression),
             legend = T, annotation_legend = T, cutree_rows = 4,  cutree_cols = 2)

a

ggsave(filename=paste0(plotdir, Sys.Date(),"_heatmap_scrna_cutcolsrows.png"),
       plot = a, device = png(), 
       width=12, 
       height=12, 
       units = "in", dpi = 300, scale = 1)

dev.off()
dev.off()

## PCA ------------
PCA_modelobject <- prcomp(scrna.df, scale=TRUE, center=TRUE)
summary(PCA_modelobject)

plotdata_cPCA <- data.frame(PCA_modelobject$x[,1:3])
plotdata_cPCA$Disease_Identity <- annotation_col$Disease_Identity

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
a <- ggplot(data = plotdata_cPCA, 
            aes(x=PC1,y=PC2, color = Disease_Identity))+
  geom_point(size=4)+ xlab(paste0("PC1 ", round(pov[1]*100, 1), "%"))+ylab(paste0("PC2 ", round(pov[2]*100,1), "%"))+
  expand_limits(x = mean(plotdata_cPCA[,"PC1"])*1.5) +theme_bw()+
  theme_journal +theme(panel.grid.major = element_line(color = "light grey",
                                                       size=0.25,
                                                       linetype = "dashed"), axis.text.x = element_text(angle = 0))+
  theme(legend.position = "right")+
  scale_color_manual(values = color, breaks = names(color))
a
ggsave(file=paste0(plotdir, as.character(Sys.Date()),"_", "PCA_scrna.png"), plot=a,
       device = "png", units = c("in"), width = 7, height = 5)


#PCA - biplot -----
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

## box plot graph ---------
colnames(cell.freq) <- case_when(colnames(cell.freq) =="B" ~ "B cells",
                                colnames(cell.freq) =="B_Plasma"  ~ "B Plasma",
                                colnames(cell.freq) =="ILC_A" ~ "ILC_A",
                                colnames(cell.freq) =="ILC_B" ~ "ILC_B",
                                colnames(cell.freq) =="Macrophage" ~ "Macrophage",
                                colnames(cell.freq) =="Macrophage_Alveolar" ~ "Macrophage_Alveolar",
                                colnames(cell.freq) =="NK" ~ "NK cells",
                                colnames(cell.freq) =="T" ~ "T helper",
                                colnames(cell.freq) =="T_Cytotoxic" ~ "T cytotoxic",
                                colnames(cell.freq) =="T_Regulatory" ~ "Tregs",
                                colnames(cell.freq) =="cDC2" ~ "cDC2",
                                colnames(cell.freq) =="cMonocyte" ~ "Classical monocytes",
                                colnames(cell.freq) =="ncMonocyte" ~ "Non-classical monocytes",
                                colnames(cell.freq) =="pDC" ~ "pDC",
                                colnames(cell.freq) =="DC_Langerhans" ~ "DC_Langerhans",
                                colnames(cell.freq) =="Mast" ~ "Mast cells",
                                colnames(cell.freq) =="cDC1" ~ "cDC1",
                                colnames(cell.freq) =="DC_Mature" ~ "DC_Mature",
                                colnames(cell.freq) =="Subject_Identity" ~ "Subject_Identity",
                                colnames(cell.freq) =="Disease_Identity" ~ "Disease_Identity")

a <- cell.freq %>% pivot_longer(where(is.numeric)) %>%
  filter(!Subject_Identity %in% c("137CO", "152CO")) %>%
  filter(!Disease_Identity %in% c("IPF")) %>%
  mutate(Disease_Identity = factor(Disease_Identity, levels = c("COPD", "Control"))) %>%
  mutate(value = log(value+1)) %>%
  ggplot(aes(value, name, fill = Disease_Identity)) +
  geom_boxplot()+ theme_bw()+
  scale_fill_manual(values = color, breaks = names(color))+
  scale_color_manual(values = color, breaks = names(color))+
  labs(y = "",      
       #y = expression("%"*CD45^"+"~cells)) +
       x= expression("%"*CD45^"+"~cells~(LOG~transformed)))+
  scale_y_discrete(limits=rev)+theme_journal+
  theme(legend.position = "top")
a
a_plot <- a+theme(legend.position = "none")
a_plot
library(ggplotify)
library(cowplot)
library(ggpubr)
legend <- as.ggplot(get_legend(a))
legend
ggsave(file=paste0(plotdir, as.character(Sys.Date()),"_", "boxplot_scrna.png"), plot=a,
       device = "png", units = c("in"), width = 7.5, height = 10.5, dpi = 300)
ggsave(file=paste0(plotdir, as.character(Sys.Date()),"_", "legend_boxplot_scrna.png"), plot=legend,
       device = "png", units = c("in"), width = 6, height = 5, dpi = 300)

# significance test of box plot -------------------------
library(ggpubr)
library(rstatix)
dt <- cell.freq %>% pivot_longer(where(is.numeric)) %>%
  mutate(Disease_Identity = factor(Disease_Identity, levels = c("COPD", "Control"))) %>%
  mutate(value = log1p(value))
my_comparisons <- list(c("COPD", "Control"))
t_test_tooth <- compare_means(value ~ Disease_Identity, group.by= "name", 
                              comparisons = my_comparisons, 
                              p.adjust.method = "fdr", 
                              method='wilcox.test', data = dt)
write_csv(t_test_tooth, file = paste0(dataoutputdir, "significance_boxplot.csv"))

# UMAP ---------------------------------------------------
library(Seurat)
a <- DimPlot(kaminski_new, group.by = "new_ID", split.by = "Disease_Identity", label = T,
             repel = T)&NoAxes()&NoLegend()&ggtitle("")
ggsave(file=paste0(plotdir, as.character(Sys.Date()),"_", "umap_scrna.png"), plot=a,
       device = "png", units = c("in"), width = 10, height = 5, dpi = 300)

#revision on UMAp ---------
#we want to create the figure 2 in paper
library(Seurat)
plotdir <- "/home/isilon/users/o_syarif/COPD machine learning/COPD_paper/Fig2/output/plots/"
datainputdir <- "/home/isilon/users/o_syarif/COPD machine learning/COPD_paper/Fig2/input/"
dataoutputdir <- "/home/isilon/users/o_syarif/COPD machine learning/COPD_paper/Fig2/output/data/"


#for publication we start from here!
#for visualisation we load the relabelled scRNA
load(paste0(datainputdir, "scrna_relabel.RData"))


library(Seurat)
color <- c( "darkgrey", "darkred")
names(color) <-  c("Control","COPD")
Idents(kaminski_new) <-"new_ID"
unique(kaminski_new$new_ID)
cols <- rep("darkgrey", 18) 
names(cols) <- unique(kaminski_new$new_ID)

sub_donor <- subset(kaminski_new, Disease_Identity=="Control")
a <- DimPlot(sub_donor, label = F,
             repel = F, cols = cols)&NoAxes()&NoLegend()&ggtitle("Control")&
  theme(plot.title = element_text(hjust=0.5));a

a_plot <-LabelClusters(
  a,
  id="ident",
  repel = TRUE,
  box = TRUE,
  fill="white",
  geom = "GeomPoint",
  position = "nearest"
);a_plot

ggsave(file=paste0(as.character(Sys.Date()),"_", "umap_scrna_donor.png"), plot=a_plot,
       device = "png", units = c("in"), width = 5, height = 5, dpi = 300)

cols <- rep("darkred", 18) 
names(cols) <- unique(kaminski_new$new_ID)
sub_copd <- subset(kaminski_new, Disease_Identity=="COPD")
a <- DimPlot(sub_copd, label = F,
             repel = F, cols = cols)&NoAxes()&NoLegend()&ggtitle("COPD")&
  theme(plot.title = element_text(hjust=0.5));a

a_plot <-LabelClusters(
  a,
  id="ident",
  repel = TRUE,
  box = TRUE,
  fill="white",
  geom = "GeomPoint",
  position = "nearest"
);a_plot

ggsave(file=paste0(as.character(Sys.Date()),"_", "umap_scrna_copd.png"), plot=a_plot,
       device = "png", units = c("in"), width = 5, height = 5, dpi = 300)

