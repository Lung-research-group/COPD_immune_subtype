# Fig 4C top


## Load the required libraries --------------
library(ggbiplot) #for ploting prcomp results with pretty circles etc
library(factoextra) #causes dplyr plyr problems !!! reload dplyr!


library(tidyverse)
library(ggpubr)
library(openxlsx)
library(ggpubr)
library(cowplot)
require(ggrepel)
library(patchwork)
library(readxl)
library(tidyr)
library(missMDA)
library(pheatmap)


## Formating ---------------
theme_set(theme_pubr(base_size=16, border = T))


man.fill <-   scale_fill_manual(values = c("Donor" = "#757b87", "COPD" = "#791812"))
man.col <-   scale_color_manual(values = c("Donor" = "#757b87", "COPD" = "#791812"))

axis <- theme(axis.line.x = element_line(color="black", size = 0.75),
              axis.line.y = element_line(color="black", size = 0.75), 
              plot.title = element_text(size = rel(1), colour="black", element_text(hjust = 0.5)),
              
              axis.title.x = element_text(margin=margin(0,10,0,0)), 
              axis.title.y = element_text(margin=margin(0,10,0,0)), 
              axis.text = element_text(colour="black"), 
              
              panel.grid = element_blank(),
              plot.caption = element_text(size = rel(0.6)))


#set the PCs to be shown
axis.x=1
axis.y=2

theme_journal2 <- theme(panel.grid.major = element_line(color = "light grey", linewidth = 0.2, linetype = "solid"), 
                        panel.grid.minor = element_line(color = "light grey", linewidth = 0.05, linetype = "solid"), 
                        plot.title = element_text(size = rel(1), colour="black", element_text(hjust = 0.5)),
                        axis.line.x = element_line(color="black", linewidth = 0.2),
                        axis.line.y = element_line(color="black", linewidth = 0.2),
                        axis.title.x = element_text(margin=margin(10,0,0,0)), 
                        axis.title.y = element_text(margin=margin(0,10,0,0)), 
                        axis.text = element_text(color="black"),
                        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.2))

symnum.args <- list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), 
                    symbols = c("***", "**", "*", "ns"))

p.adj <- function(x) {p.adjust(x, method = "BH")}

signif  <- function(x) {ifelse(x > -log10(0.05), "yes", "no")}



PCA3 <- function(dataframe, title) {
  
  p <- dataframe %>% column_to_rownames(var = "Unique_Sample_ID") %>% select(where(is.numeric))
  #x <<- colSums(is.na(p))< nrow(p)*0.3 #remove columns with >30% missing
  #q <- p[,x]
  #y <<- rowSums(is.na(p))< ncol(p)*0.3 #remove rows with >30% missing
  q <- p
  
  q <- q %>% mutate(across(.cols = everything(), ~ ifelse(is.infinite(.x), 0, .x)))
  
  q <- imputePCA(q, 5)$completeObs #imput the missing datasets
  
  df1 <<- prcomp(q, scale = T, center = T) 
  pcVar <- summary(df1)
  biplot(df1, scale = 0)
  rota <<- df1$rotation
  varPC1 <<- round(pcVar$importance[2,1], digits = 3)
  varPC2 <<- round(pcVar$importance[2,2], digits = 3)
  varPC3 <<- round(pcVar$importance[2,3], digits = 3)
  
  PC <<- df1$x %>% as.data.frame() %>% rownames_to_column(var = "Unique_Sample_ID") %>% left_join(select(dataframe, !where(is.numeric)))
  
  p1 <<- ggplot(PC,aes(x=PC1,y=PC2,col=COPD_subclass))+
    geom_point(size=4,alpha=1)+
    scale_color_manual(name = "Subtype", labels=c("A" = 'Adaptive Immune Severe', "B" = 'Immune Mild'), values=c("#661510", "#e34e45"))+
    labs(x = paste0("PC1 (", varPC1*100, "%)"), y = paste0("PC2 (", varPC2*100, "%)"), title = title)+
    theme(legend.position = "none")+
    theme(panel.grid.major = element_line(size = 0.75))+ theme_journal
  
  p2 <<- fviz_pca_biplot(df1, #repel = TRUE, 
                         axes = c(1,2), 
                         habillage= PC$COPD_subclass, #label = "var",
                         palette = c("#661510", "#e34e45"),
                         mean.point = FALSE, pointsize = 4, pointshape = 19, 
                         label ="var", select.var = list(contrib = 10), 
                         labelsize = 4, arrowsize = 0.1, col.var = "#252525", repel = TRUE,
                         title = paste0(title, " PCA - Biplot")
  )+theme(legend.position = "right") + theme_journal
  
  p1+p2
  
}

PCA2 <- function(dataframe, title) {
  
  p <- dataframe %>% column_to_rownames(var = "Unique_Sample_ID") %>% select(where(is.numeric))
  
  x <<- colSums(is.na(p))< nrow(p)*0.3 #remove columns with >30% missing
  q <- p[,x]
  y <<- rowSums(is.na(p))< ncol(p)*0.3 #remove rows with >30% missing
  q <- q[y,]
  #q <- log10(q) 
  
  q <- q %>% mutate(across(.cols = everything(), ~ ifelse(is.infinite(.x), 0, .x)))
  
  q <- imputePCA(q, 5)$completeObs #imput the missing datasets
  
  df1 <<- prcomp(q, scale = T, center = T) 
  pcVar <- summary(df1)
  biplot(df1, scale = 0)
  rota <<- df1$rotation
  varPC1 <<- round(pcVar$importance[2,1], digits = 3)
  varPC2 <<- round(pcVar$importance[2,2], digits = 3)
  varPC3 <<- round(pcVar$importance[2,3], digits = 3)
  
  PC <<- df1$x %>% as.data.frame() %>% rownames_to_column(var = "Unique_Sample_ID") %>% left_join(select(dataframe, !where(is.numeric)))
  
  p1 <<- ggplot(PC,aes(x=PC1,y=PC2,col=COPD_subclass))+
    geom_point(size=4,alpha=1)+
    #geom_label_repel()+
    scale_color_manual(name = "Subtype", labels=c("A" = 'Adaptive Immune Severe', "B" = 'Immune Mild'), values=c("#661510", "#e34e45"))+
    labs(x = paste0("PC1 (", varPC1*100, "%)"), y = paste0("PC2 (", varPC2*100, "%)"), title = title)+
    theme(legend.position = "none")+
    theme(panel.grid.major = element_line(size = 0.75))+ theme_journal
  
  p2 <<- fviz_pca_biplot(df1, #repel = TRUE, 
                         axes = c(1,2), 
                         habillage= PC$COPD_subclass, #label = "var",
                         palette = c("#757b87", "#791812"),
                         mean.point = FALSE, pointsize = 4, pointshape = 19, 
                         label ="var", select.var = list(contrib = 10), 
                         labelsize = 4, arrowsize = 0.1, col.var = "#252525", repel = TRUE,
                         title = paste0(title, " PCA - Biplot")
  )+theme(legend.position = "right") + theme_journal
  
  p1+p2
  
}



# Import data -------------------
Parameter_labels <- read_excel("iScience_input_data/Parameter_labels.xlsx", 1)
parameter_labeller <- Parameter_labels$expression
names(parameter_labeller) <- Parameter_labels$Parameter_Name

facs.per.trans <- read.xlsx("iScience_input_data/hFACS_Donor_COPD_per_LOG_x1c_full.xlsx") %>% mutate(Unique_Sample_ID = as.factor(Unique_Sample_ID))
elisa.serum.25.c1 <- read.xlsx("iScience_input_data/hELISAS_Donor_COPD_conc_LOG_0quarterofmin_cohort1.xlsx") %>% mutate(Unique_Sample_ID = as.factor(Unique_Sample_ID))

COPD_subclass <- read.xlsx("iScience_input_data/COPD_subclass.xlsx") %>% 
  mutate(COPD_subclass = as.factor(COPD_subclass)) %>% mutate(Unique_Sample_ID = as.factor(Unique_Sample_ID))

double.d <- c("Lymphocytes", # parental pop
              "Macrophages", # parental pop
              #"Monocytes", # parental pop
              "PMNL", # parental pop
              "Mast", # parental pop
              "Macs_CD14med", "Macs_CD14hi", # parental pop
              "DC_CD209neg_CD11cpos", 
              "DC_CD209pos_CD11cneg", 
              "DC_CD209pos_CD11cpos") 



# Fig 4C -----------
## FACS Subtypes + plasma -----------
temp <- elisa.serum.25.c1 %>% mutate(across(where(is.numeric), log10)) 
data <- COPD_subclass %>% left_join(facs.per.trans)%>%
  select(!double.d) %>%
  left_join(temp) 

#View(PC)
PCA2(data, "%CD45 cells + plasma")


fviz_pca_biplot(df1, title = NULL, 
                axes = c(axis.x, axis.y), 
                habillage= PC$COPD_subclass, #label = "var",
                palette = c("#661510", "#e34e45"),
                mean.point = FALSE, pointsize = 4, pointshape = 19, 
                label ="var", select.var = list(contrib = 14), #addEllipses=TRUE, ellipse.level=0.95
                labelsize = 4, arrowsize = 0.1, col.var = "#252525", repel = TRUE,
                subtitle = expression("Percentage CD45"^"+"~"cells and plasma")
)+
  expand_limits(x = c(min(df1$x[,paste0("PC", axis.x)])*1.1, max(df1$x[,axis.x])*1.1),
                y = c(min(df1$x[,paste0("PC", axis.y)])*1.1, max(df1$x[,axis.y])*1.1)) +  
  labs(x = paste0("PC", axis.x, " (", round(summary(df1)$importance[2,axis.x], digits = 3)*100, "%)"),
       y = paste0("PC", axis.y, " (", round(summary(df1)$importance[2,axis.y], digits = 3)*100, "%)"))+
  theme_pubr(base_size = 16)+
  
  theme(legend.position = "none") + theme_journal2

save_plot(paste0("plots/", Sys.Date(), "_COPD_4C", ".png"), 
          plot = ggplot2::last_plot(), base_height = 7, base_width = 7, dpi = 300)
save_plot(paste0("plots/", Sys.Date(), "_COPD_4C", ".pdf"), 
          plot = ggplot2::last_plot(), base_height = 7, base_width = 7, dpi = 300)

ggsave(paste0(output_folder, model_name,"__cPCA_PC12_biplot.png"),plot = Biplot_cPCA_PC12 , device = png(), width=w+5, 
       height=h+7, units = "mm", dpi = 900); dev.off()
ggsave(paste0(output_folder, model_name,"__cPCA_PC12_biplot.pdf"),plot = Biplot_cPCA_PC12 , device = pdf(), width=w+5, 
       height=h+7, units = "mm", dpi = 900); dev.off()

# Fig 4B -----------
## FACS Subtypes + plasma -----------
temp <- elisa.serum.25.c1 %>% mutate(across(where(is.numeric), log10)) 
data <- COPD_subclass %>% left_join(facs.per.trans)%>%
  select(!double.d) #%>%
#left_join(temp) 

#View(PC)
PCA2(data, "%CD45 cells")


fviz_pca_biplot(df1, title = NULL, 
                axes = c(axis.x,axis.y), 
                habillage= PC$COPD_subclass, #label = "var",
                palette = c("#661510", "#e34e45"),
                mean.point = FALSE, pointsize = 4, pointshape = 19, 
                label ="var", select.var = list(contrib = 10), #addEllipses=TRUE, ellipse.level=0.95
                labelsize = 4, arrowsize = 0.1, col.var = "#252525", repel = TRUE,
                subtitle = expression("Percentage CD45"^"+"~"cells and plasma")
)+
  expand_limits(x = c(min(df1$x[,paste0("PC", axis.x)])*1.1, max(df1$x[,axis.x])*1.1),
                y = c(min(df1$x[,paste0("PC", axis.y)])*1.1, max(df1$x[,axis.y])*1.1)) +  
  labs(x = paste0("PC", axis.x, " (", round(summary(df1)$importance[2,axis.x], digits = 3)*100, "%)"),
       y = paste0("PC", axis.y, " (", round(summary(df1)$importance[2,axis.y], digits = 3)*100, "%)"))+
  theme_pubr(base_size = 16)+
  
  theme(legend.position = "none") + theme_journal2

save_plot(paste0("plots/", Sys.Date(), "_COPD_4B", ".png"), 
          plot = ggplot2::last_plot(), base_height = 7, base_width = 7, dpi = 300)
save_plot(paste0("plots/", Sys.Date(), "_COPD_4B", ".pdf"), 
          plot = ggplot2::last_plot(), base_height = 7, base_width = 7, dpi = 300)
