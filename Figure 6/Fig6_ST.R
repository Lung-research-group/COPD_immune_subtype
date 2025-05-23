#we will prepare for cibersortx
library(dplyr)
library(reshape2)
library(Seurat)
library(stringr)
library(ggplot2)

#process the cibersortx results ------------------
dir <- paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/")
cibersort <- read.csv(file =paste0(dir, "input data/CIBERSORTx_Job22_sub50_Control_COPD_15pops_parenchyma_sub.csv"))

#filtering out non significant deconvolutions
pop <- cibersort %>%
  filter(P.value <0.05)

names_pop <- colnames(pop)[-which(colnames(pop) %in% c("P.value", "Correlation","RMSE", "Absolute.score..sig.score."))]
names_pop

library(tibble)
pop <- column_to_rownames(pop, var="Mixture")
pop <- pop[-which(colnames(pop) %in% c("P.value", "Correlation","RMSE", "Absolute.score..sig.score."))]

#the cell types of interest
immune <- c("Classical.monocytes", "Non.classical.monocytes", "Macrophage","Macrophage_Alveolar",
            "Mast.cells", "pDC","cDC2","cDC1", "DC_Mature", 
            "B.cells", "B.Plasma", "T.cytotoxic"  , 
            "Tregs","T.helper",               
            "NK.cells")

#descriptive of the data
melt <- melt(pop)
melt <- na.omit(melt)
melt$variable <- factor(melt$variable, levels = immune)
a <- ggplot(data=melt, aes(x=value))+
  geom_histogram()+facet_wrap(~variable, scales = "fixed")

a

#perform transformation
melt$transform <- log10(melt$value+0.001)

#check how transformation change the data distribution
a <- ggplot(data=melt, aes(x=transform))+
  geom_histogram()+facet_wrap(~variable, scales = "fixed")

a

#perform transformation on original data
library(tidyverse)
dt <- as.data.frame(sapply(pop, function(num) log10(num+0.001)))
rownames(dt) <- rownames(pop)

# patient level data  -----------------------------------------
data <- read.csv(file =paste0(dir, "input data/metadata_nanostring.csv"))

#combine the metadata from original nanostring
patient <- cbind.data.frame("patient"=data[["patient_ID"]],
                            "location"=data[["location"]],
                            "ROI"=data[["ROI"]]) # addlocation
rownames(patient) <- data$ROI
patient <- patient[which(rownames(patient) %in% rownames(pop)),]

#continue
pati <- cbind.data.frame(dt,
                         patient)
melt <- melt(pati)
melt$patient_variable <- paste0(melt$patient, "_", melt$variable)
melt <- melt%>% group_by(patient_variable) %>%
  summarise(mean = mean(value))
temp <- as.data.frame(str_split_fixed(melt$patient_variable, "_", 2))
melt$patient <- temp$V1
melt$variable <- temp$V2

meta <- cbind.data.frame("gold_broad"=data[["gold_broad"]], 
                         "patient"=data[["patient_ID"]], 
                         "emphysema"=data[["emphysema_whole"]])

meta <-unique(meta)

merge <- merge(meta, melt, by="patient")

library(tidyr)
df <- pivot_wider(merge[, c("patient", "mean", "variable", "gold_broad","emphysema")],
                  names_from = "variable", values_from = "mean")
df$X <- NULL #remove additional column

#relabel the cells
colnames(df) <- case_when(colnames(df) == "patient"~ "patient",
                          colnames(df) == "gold_broad"~ "gold_broad",
                          colnames(df) == "emphysema"~ "emphysema",
                          colnames(df) == "T.cytotoxic"~ "T cytotoxic",
                          colnames(df) == "Classical.monocytes"~ "Classical monocytes",
                          colnames(df) == "pDC"~ "pDC",
                          colnames(df) == "T.helper"~ "T helper",
                          colnames(df) == "DC_Mature"~ "DC_Mature",
                          colnames(df) == "B.Plasma"~ "B Plasma",
                          colnames(df) == "B.cells"~ "B cells",
                          colnames(df) == "Macrophage_Alveolar"~ "Macrophage_Alveolar",
                          colnames(df) == "Non.classical.monocytes"~ "Non classical monocytes",
                          colnames(df) == "Macrophage"~ "Macrophage",
                          colnames(df) == "cDC2"~ "cDC2",
                          colnames(df) == "Tregs"~ "Tregs",
                          colnames(df) == "Mast.cells"~ "Mast cells",
                          colnames(df) == "cDC1"~ "cDC1",
                          colnames(df) == "NK.cells"~ "NK cells")


# kmeans clustering ----------------------
# Compute k-means with k = 2
set.seed(123)
rownames(df) <- df$patient
dfx <- scale(df[which(df$gold_broad %in% c("GOLD1-2", "GOLD3-4")), c(4:18)]) # only to take the GOLD 3-4
km.res <- kmeans(dfx, 2, nstart = 50)

# Print the results
print(km.res)

dd <- cbind(df[which(df$gold_broad %in% c("GOLD1-2","GOLD3-4")),], cluster = km.res$cluster)
head(dd)
rownames(dd) <- dd$patient
table(dd$gold_broad, dd$cluster)
cluster2 <- dd$patient[which(dd[, "cluster"] %in% c(1))] # we swap the naming to align with subgroup FACS results
cluster1 <- dd$patient[which(dd[, "cluster"] %in% c(2))]
include <- c(cluster1, cluster2)

## bi-plot PCA (fig5B) -----------------------------
PCA_modelobject <- prcomp(df[,c(4:18)], scale=T, center=T)
summary(PCA_modelobject)

plotdata_cPCA <- data.frame(PCA_modelobject$x[,1:3])

#we take patients ID from MIC.R (available in PC)
GOLD12 <- unique(df$patient[which(df$gold_broad %in% c("GOLD1-2"))])
plotdata_cPCA$patient <- df$patient[which(df$gold_broad %in% c("GOLD1-2","GOLD3-4"))]
plotdata_cPCA$subgroup <- ifelse(plotdata_cPCA$patient %in% cluster2, "cluster_2","cluster_1") #clust 1= severe clust 2 mild
plotdata_cPCA$gold_broad <- ifelse(plotdata_cPCA$patient %in% GOLD12, "GOLD1-2","GOLD3-4")
pov <- PCA_modelobject$sdev^2/sum(PCA_modelobject$sdev^2)
barplot(pov)
theme_journal <- theme(axis.line.x = element_line(color="black", size = 0.5),
                       axis.line.y = element_line(color="black", size = 0.5))+
  theme(axis.text.x = element_text(colour = "black", size = 16, angle = 0, hjust=1, vjust=1))+
  theme(axis.text.y = element_text(colour = "black", size = 16))+
  theme(axis.title = element_text(size = 16, colour = "black"))+
  theme(legend.position="bottom")+
  theme(text = element_text(size = 16))+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16))

plotdata_cPCA$subgroup <- factor(plotdata_cPCA$subgroup, levels = c("cluster_1","cluster_2"))
plotdata_cPCA$subgroup <- ifelse(plotdata_cPCA$subgroup == "cluster_1", "cluster 1", "cluster 2")
color <- c( "darkred", "salmon")
names(color) <- c("cluster 1","cluster 2")
cols <- c("#CD5C5C","#791812")
names(cols) <- c("GOLD1-2","GOLD3-4")
library(RColorBrewer)

library(ggplot2)
library(ggpubr)
theme_set(theme_pubr(base_size=6, border = T, legend = "none"))
theme_journal <- theme(panel.grid.major = element_line(color = "light grey", linewidth = 0.2, linetype = "solid"),
                       panel.grid.minor = element_line(color = "light grey", linewidth = 0.05, linetype = "solid"),
                       plot.title = element_text(size = rel(1)),
                       axis.line.x = element_line(color="black", linewidth = 0.2),
                       axis.line.y = element_line(color="black", linewidth = 0.2),
                       axis.text = element_text(color="black", size = 6),
                       panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.2)
)
library(factoextra)
Biplot_cPCA <- fviz_pca_biplot(PCA_modelobject, title = "PCA-Biplot", 
                               habillage= plotdata_cPCA$subgroup,
                               mean.point = FALSE, pointsize = 1, pointshape = 19,
                               label ="var", select.var = list(contrib = 10), #addEllipses=TRUE, ellipse.level=0.95
                               col.var = "black",labelsize = 2, arrowsize = 0.12, repel = TRUE) +  #can't color by contrib since uses already color_scale for PCA scores
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


ggsave(paste0(plotdir, "ST__cPCA_biplot_PCA_final2.png"),
       plot = Biplot_cPCA , 
       device = png(), width=65, height=70, units = "mm", dpi = 900)
dev.off()

## violin plots (fig5D) =====
nnaostring_metadata <- read.csv(file = paste0(dir, "input data/metadata_clinical_data.csv"))
nnaostring_metadata$X <- NULL
meta <- cbind.data.frame("LAA950_whole"=data[["LAA950_whole"]],
                         "LAA950_lobe"=data[["LAA950_lobe"]],
                         "patient"=data[["patient_ID"]])

patient_ID <- unique(meta$patient)
overlap <- intersect(nnaostring_metadata$patient, meta$patient)
nnaostring_metadata <- nnaostring_metadata[which(nnaostring_metadata$patient %in% patient_ID),]
nnaostring_metadata <- nnaostring_metadata[, c("patient", "FEV1.pred", 
                                               "FEV1.FVC" ,"Pack.years")] # "progression.type"
merge <- merge(meta, nnaostring_metadata, by="patient")
merge <-unique(merge)
meta <- merge
merge <- merge(meta, df, by="patient")

gold_12 <- merge[which(merge$gold_broad %in% c("GOLD1-2", "GOLD3-4")),-which(colnames(merge) %in% c("LAA950_lobe"))]
gold_12 <- gold_12[which(gold_12$patient %in% include),]
gold_12 <- unique(gold_12)
gold_12$status <- ifelse(gold_12$patient %in% cluster2, "cluster_2","cluster_1")

library(ggpubr)
library(rstatix)
melt <- melt(gold_12)
my_comparisons <- list(c("cluster_1","cluster_2"))
t_test_tooth <- compare_means(value ~ status, group.by= "variable", 
                              comparisons = my_comparisons, 
                              p.adjust.method = "fdr", 
                              method='wilcox.test', data = melt)

zeze <- melt %>%group_by(variable) %>%
  mutate(position = max(value, na.rm = T))
zeze <- zeze[c("variable", "position")]
zeze <- unique(zeze)
#t_test_tooth$y.position <- NULL
t_test_tooth$position <- NULL
t_test_tooth <- merge(t_test_tooth, zeze, by=c("variable"))
t_test_tooth$cut <- case_when(t_test_tooth$p.adj >0.05 ~ "ns",
                              t_test_tooth$p.adj <=0.05 & t_test_tooth$p.adj >0.01 ~ "*",
                              t_test_tooth$p.adj <=0.01 & t_test_tooth$p.adj >0.001 ~ "**",
                              t_test_tooth$p.adj <=0.001 & t_test_tooth$p.adj > 0.0001~ "***",
                              t_test_tooth$p.adj <=0.0001 ~ "****")
t_test_tooth$stat <- paste0(t_test_tooth$group1, "_", t_test_tooth$group2)
t_test_tooth$y.position <- t_test_tooth$position +5 #0.3 or 5

uni <- t_test_tooth[, c("variable", "position", "group1")]

# Now filter out the cell types that we want
pop <- c("Mast cells", "Macrophage", "pDC", 
         "T helper")
melt <- filter(melt, variable %in% pop)
uni <- filter(uni, variable %in% pop)
t_test_tooth_01 <- filter(t_test_tooth, variable %in% pop)


library(ggbeeswarm)
library(scales)


cols <- c("#CD5C5C","#791812")
names(cols) <- c("GOLD1-2","GOLD3-4")

color <- c( "darkred", "salmon")
names(color) <- c("cluster_1","cluster_2")

color <- c( "darkred", "salmon")
names(color) <- c("cluster 1","cluster 2")

#for publication
pop <- c("Mast cells", "Macrophage", "pDC", 
         "T helper")

symnum.args <- list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), symbols = c("***", "**", "*", "ns"))

melt$status <- ifelse(melt$status == "cluster_1", "cluster 1", "cluster 2")
labelss <- c( "cluster\n1", "cluster\n2")

p <-ggplot(melt, aes(x=status, y=value)) +
  geom_violin(color="white", alpha=0.2, lwd= 1, aes(group=status, fill=status))+
  geom_dotplot(dotsize = 2, binaxis="y", stackdir="center", 
               aes(fill=status)) +
  stat_compare_means(aes(label = paste0(after_stat(p.signif))), size = 8,  
                     symnum.args = symnum.args, label.x.npc = "center",
                     label.y.npc = "top")+
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.25)))+
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))+
  facet_wrap(~variable, scales = "free_y", ncol = 4)+
  scale_fill_manual(values=color, breaks = names(color))+
  scale_color_manual(values=color, breaks = names(color))+ xlab("")+
  scale_x_discrete(label=labelss)+ylab("abundance score\n(LOG transformed)")+
  theme(strip.text.x = element_text(size = 16),
        axis.text = element_text(size=16),
        axis.title.y=element_text(size=16))
p
p <- p+theme(legend.position = "none")
p

library(cowplot)
save_plot(filename=paste0(plotdir, Sys.Date(),"rev01_postave_COPD_Control_GOLDall_ p value_individual analysis_cibersortx_parenchyma_publication",".png"),
          plot =p , base_height = 5.25, base_width = 12.5, dpi = 600)


## bar chart (fig5E) -------------------
nnaostring_metadata <- read.csv(file = paste0(dir, "input data/metadata_clinical_data.csv"))
nnaostring_metadata$X <- NULL
colnames(nnaostring_metadata)
meta <- cbind.data.frame("LAA950_whole"=data[["LAA950_whole"]],
                         "LAA950_lobe"=data[["LAA950_lobe"]],
                         "patient"=data[["patient_ID"]])

patient_ID <- unique(meta$patient)
overlap <- intersect(nnaostring_metadata$patient, meta$patient)
nnaostring_metadata <- nnaostring_metadata[which(nnaostring_metadata$patient %in% patient_ID),]
nnaostring_metadata <- nnaostring_metadata[, c("patient", "FEV1.pred", 
                                               "FEV1.FVC" ,"Pack.years")] # "progression.type"
merge <- merge(meta, nnaostring_metadata, by="patient")
merge <-unique(merge)
meta <- merge
merge <- merge(meta, df, by="patient")

gold_12 <- merge[which(merge$gold_broad %in% c("GOLD1-2", "GOLD3-4")),-which(colnames(merge) %in% c("LAA950_lobe"))]
gold_12 <- gold_12[which(gold_12$patient %in% include),]
gold_12 <- unique(gold_12)
gold_12$status <- ifelse(gold_12$patient %in% cluster2, "cluster_2","cluster_1")
gold_12$status_gap <-ifelse(gold_12$gold_broad %in% c("GOLD1-2")& gold_12$LAA950_whole >10, "GOLD1-2 progress", "rest")

library(ggpubr)
library(rstatix)
melt <- melt(gold_12)
lala <- melt[, c("gold_broad", "status")]

lala <- as.data.frame(table(lala))
lala$name <- "Gold Broad Category"
cols <-c("#D9B926", "#D02F3B")
names(cols) <- c("GOLD1-2","GOLD3-4")

p2 <- ggplot(lala, aes(fill=gold_broad, y=Freq, x=status)) + 
  geom_bar(position="fill", stat="identity")+
  scale_y_continuous(labels = scales::percent) +
  scale_x_discrete(labels=c("cluster\n1", "cluster\n2"))+
  labs(y = "Relative proportion", x = "")+
  scale_fill_manual(values=cols, breaks = names(cols))+
  facet_wrap(~name, scales = "free_y", ncol = 1)+
  theme(strip.text.x = element_text(size = 16),
      axis.text = element_text(size=16),
      axis.title.y=element_text(size=16),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.key.spacing.x = unit(0.2, "in") )
p2 <- p2+theme(legend.position = "top") #set to 'none' for plotting only bar chart
p2
save_plot(filename=paste0(plotdir, Sys.Date(),"rev01_postave_COPD_Control_GOLDall_ p value_individual analysis_cibersortx_parenchyma_publication_cells_byclusters_FEV1_dotplot.png"),
          plot = p2, base_height = 4.5, base_width = 3.5, dpi = 600)
legend <- as.ggplot(get_legend(p2))
legend
ggsave(filename = paste0(plotdir, Sys.Date(),"legend_rev01_postave_COPD_Control_GOLDall_ p value_individual analysis_cibersortx_parenchyma_publication_cells.png"),
       plot=legend,
       device = "png", width = 3, 
       height = 1, dpi = 600,
       units = c("in"), scale = 1)


## violin plot (fig5F) -----------------------------
nnaostring_metadata <- read.csv(file = paste0(dir, "input data/metadata_clinical_data.csv"))
nnaostring_metadata$X <- NULL
meta <- cbind.data.frame("LAA950_whole"=data[["LAA950_whole"]],
                         "LAA950_lobe"=data[["LAA950_lobe"]],
                         "patient"=data[["patient_ID"]])

patient_ID <- unique(meta$patient)
overlap <- intersect(nnaostring_metadata$patient, meta$patient)
nnaostring_metadata <- nnaostring_metadata[which(nnaostring_metadata$patient %in% patient_ID),]
nnaostring_metadata <- nnaostring_metadata[, c("patient", "FEV1.pred", 
                                               "FEV1.FVC" ,"Pack.years")] # "progression.type"
merge <- merge(meta, nnaostring_metadata, by="patient")
merge <-unique(merge)
meta <- merge
merge <- merge(meta, df, by="patient")

gold_12 <- merge[which(merge$gold_broad %in% c("GOLD1-2", "GOLD3-4")),-which(colnames(merge) %in% c("LAA950_lobe"))]
gold_12 <- gold_12[which(gold_12$patient %in% include),]
gold_12 <- unique(gold_12)
gold_12$status <- ifelse(gold_12$patient %in% cluster2, "cluster_2","cluster_1")
gold_12$status_gap <-ifelse(gold_12$gold_broad %in% c("GOLD1-2")& gold_12$LAA950_whole >10, "GOLD1-2 progress", "rest")

library(ggpubr)
library(rstatix)
melt <- melt(gold_12)
my_comparisons <- list(c("cluster_1","cluster_2"))
t_test_tooth <- compare_means(value ~ status, group.by= "variable", 
                              comparisons = my_comparisons, 
                              p.adjust.method = "fdr", 
                              method='wilcox.test', data = melt)

zeze <- melt %>%group_by(variable) %>%
  mutate(position = max(value, na.rm = T))
zeze <- zeze[c("variable", "position")]
zeze <- unique(zeze)

t_test_tooth$position <- NULL
t_test_tooth <- merge(t_test_tooth, zeze, by=c("variable"))
t_test_tooth$cut <- case_when(t_test_tooth$p.adj >0.05 ~ "ns",
                              t_test_tooth$p.adj <=0.05 & t_test_tooth$p.adj >0.01 ~ "*",
                              t_test_tooth$p.adj <=0.01 & t_test_tooth$p.adj >0.001 ~ "**",
                              t_test_tooth$p.adj <=0.001 & t_test_tooth$p.adj > 0.0001~ "***",
                              t_test_tooth$p.adj <=0.0001 ~ "****")
t_test_tooth$stat <- paste0(t_test_tooth$group1, "_", t_test_tooth$group2)
t_test_tooth$y.position <- t_test_tooth$position +5 #0.3 or 5

uni <- t_test_tooth[, c("variable", "position", "group1")]

# Now filter LAA950 
melt <- filter(melt, variable == "LAA950_whole")
uni <- filter(uni, variable == "LAA950_whole")
t_test_tooth_01 <- filter(t_test_tooth, variable == "LAA950_whole")

#visualisation 
melt$status <- ifelse(melt$status == "cluster_1", "cluster 1", "cluster 2")
color <- c( "darkred", "salmon")
names(color) <- c("cluster 1","cluster 2")
labelss <- c( "cluster\n1", "cluster\n2")

#colored by clusters LAA950
p2 <-ggplot(melt, aes(x=status, y=value, group=status)) +
  geom_violin(color="white",alpha=0.2, lwd= 1, aes(fill=status))+
  geom_dotplot(dotsize = 2, binaxis="y", stackdir="center",stackgroups = T,binpositions = "all",
               aes(x=status, y=value,
                   fill=status), inherit.aes = F) +
  stat_compare_means(aes(label = paste0(after_stat(p.signif))), size = 8,  
                     symnum.args = symnum.args, label.x.npc = "center",
                     label.y.npc = "top")+
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.25)))+
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))+
  facet_wrap(~variable, scales = "free_y", ncol = 4)+
  scale_fill_manual(values=color, breaks = names(color))+
  scale_color_manual(values=color, breaks = names(color))+
  scale_x_discrete(label=labelss)+ylab("LAA950 whole lung")+xlab("")+
  theme(strip.text.x = element_text(size = 16),
        axis.text = element_text(size=16),
        axis.title.y=element_text(size=16))
p2
p_plot <-p2+theme(legend.position="none")
save_plot(filename=paste0(plotdir, Sys.Date(),"rev01_postave_COPD_Control_GOLDall_ p value_individual analysis_cibersortx_parenchyma_publication_cells_byclusters.png"),
          plot = p_plot, base_height = 5.25, base_width = 3.5, dpi = 600)

## correlation plots (fig5G) ====
meta <- cbind.data.frame("LAA950_whole"=data[["LAA950_whole"]],
                         "LAA950_lobe"=data[["LAA950_lobe"]],
                         "patient"=data[["patient_ID"]])

meta <-unique(meta)
patient_ID <- unique(meta$patient)

overlap <- intersect(nnaostring_metadata$patient, meta$patient)
nnaostring_metadata <- nnaostring_metadata[which(nnaostring_metadata$patient %in% patient_ID),]
nnaostring_metadata <- nnaostring_metadata[, c("patient","Pack.years" )] #"FEV1.pred", "FEV1.FVC", "Pack.years"
meta <- merge(meta, nnaostring_metadata, by="patient")
meta <-unique(meta)

melt <- melt(df)
merge <- left_join(melt, meta, by="patient")

#GOLD1-2/3-4
gold_12 <- merge[which(merge$gold_broad %in% c("GOLD1-2","GOLD3-4")),-which(colnames(merge) %in% c("LAA950_lobe", "LAA950_whole"))]
gold_12$status <- ifelse(gold_12$patient %in% cluster2,  "cluster_2","cluster_1")

library(ggpubr)

cols <-c("#D9B926", "#D02F3B")
names(cols) <- c("GOLD1-2","GOLD3-4")

immune <- c("Classical monocytes", "Non classical monocytes", "Macrophage","Macrophage_Alveolar",
            "Mast cells", "pDC","cDC2","cDC1", "DC_Mature", 
            "B cells", "B Plasma", "T cytotoxic"  , 
            "Tregs","T helper",               
            "NK cells")

merge$variable <- factor(merge$variable, levels = immune)
merge$LAA950_whole_t <- log10(merge$LAA950_whole)

sp <- ggplot(merge[which(merge$gold_broad %in% c("GOLD1-2", "GOLD3-4") &
                           merge$variable %in%
                           c("Mast cells","Macrophage",
                             "pDC","Classical monocytes")),], aes(LAA950_whole_t, value)) + 
  geom_point(aes(color=gold_broad), size=6) + 
  facet_wrap(~variable, scales = "free", ncol = 4)+
  theme(strip.text = element_text(size = 18),
        axis.text.x = element_text(size = 14, angle = 0, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        text = element_text(size = 14))+ggtitle("")+ylab("abundance score\n(LOG transformed)")+xlab("LAA950 whole lung (LOG transformed)")+
  scale_color_manual(values = cols, breaks = names(cols))+
  geom_smooth(method = "lm")+theme_bw(base_size = 16)+
  theme(axis.text= element_text(hjust=0.5, size = 18),
        axis.title= element_text(size = 20), strip.text = element_text(
          size = 20), legend.title = element_blank(),legend.key.spacing.y = unit(0.1, "in"))+ggtitle("")
sp

sp2 <-sp+stat_cor(method = "spearman", cor.coef.name="R", size=6,
                 label.y.npc="top", label.x.npc = "left")

sp2
sp_plot <- sp+theme(legend.position = "none")
sp_plot

library(ggplotify)
library(cowplot)
legend <- as.ggplot(get_legend(sp))
legend

ggsave(file = paste0(plotdir, as.character(Sys.Date()),"_", "corr_rev01_postave_COPD_Control_correlation plots_nanostring_median_cibersortx_publication.png"),
       plot = sp2,scale=1,
       width=20, height=6,dpi = 600)

ggsave(file = paste0(plotdir, as.character(Sys.Date()),"_", "legend_rev01_postave_COPD_Control_correlation plots_nanostring_median_cibersortx_publication.png"),
       plot = legend,scale=1,
       width=1.5, height=1,dpi = 600)

#end------------------------------------
