# Fig 2 + S5 & S6 Cytokines ----------------

# Setup ----------------
dir.create("Plots")

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

PCA2 <- function(dataframe, title) {
  
  p <- dataframe %>% column_to_rownames(var = "Unique_Sample_ID") %>% select(where(is.numeric))
  
  x <<- colSums(is.na(p))< nrow(p)*0.3 #remove columns with >30% missing
  q <- p[,x]
  y <<- rowSums(is.na(p))< ncol(p)*0.3 #remove rows with >30% missing
  q <- q[y,]
  q <- log10(q) 
  
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
  
  p1 <<- ggplot(PC,aes(x=PC1,y=PC2,col=Diagnosis))+
    geom_point(size=4,alpha=1)+
    man.col+
    labs(x = paste0("PC1 (", varPC1*100, "%)"), y = paste0("PC2 (", varPC2*100, "%)"), title = title)+
    theme(legend.position = "none")+
    theme(panel.grid.major = element_line(size = 0.75))+ theme_journal2
  
  p2 <<- fviz_pca_biplot(df1, #repel = TRUE, 
                         axes = c(1,2), 
                         habillage= PC$Diagnosis, #label = "var",
                         palette = c("#757b87", "#791812"),
                         mean.point = FALSE, pointsize = 4, pointshape = 19, 
                         label ="var", select.var = list(contrib = 10), 
                         labelsize = 4, arrowsize = 0.1, col.var = "#252525", repel = TRUE,
                         title = paste0(title, " PCA - Biplot")
  )+theme(legend.position = "right") + theme_journal2
  
  p1+p2
  
}

# Import data -------------------
Parameter_labels <- read_excel("iScience_input_data/Parameter_labels.xlsx", 1)
parameter_labeller <- Parameter_labels$expression
names(parameter_labeller) <- Parameter_labels$Parameter_Name

elisa.lung.25.full <- read.xlsx("iScience_input_data/hELISAL_Donor_COPD_conc_tprotNorm_OEME_LOG_0quarterofmin_full.xlsx") %>% mutate(Unique_Sample_ID = as.factor(Unique_Sample_ID))
elisa.lung.25.c1 <- read.xlsx("iScience_input_data/hELISAL_Donor_COPD_conc_tprotNorm_OEME_LOG_0quarterofmin_cohort1.xlsx") %>% mutate(Unique_Sample_ID = as.factor(Unique_Sample_ID))
elisa.lung.25.c2 <- read.xlsx("iScience_input_data/hELISAL_Donor_COPD_conc_tprotNorm_OEME_LOG_0quarterofmin_cohort2.xlsx") %>% mutate(Unique_Sample_ID = as.factor(Unique_Sample_ID))

elisa.serum.25.full <- read.xlsx("iScience_input_data/hELISAS_Donor_COPD_conc_LOG_0quarterofmin_full.xlsx") %>% mutate(Unique_Sample_ID = as.factor(Unique_Sample_ID))
elisa.serum.25.c1 <- read.xlsx("iScience_input_data/hELISAS_Donor_COPD_conc_LOG_0quarterofmin_cohort1.xlsx") %>% mutate(Unique_Sample_ID = as.factor(Unique_Sample_ID))
elisa.serum.25.c2 <- read.xlsx("iScience_input_data/hELISAS_Donor_COPD_conc_LOG_0quarterofmin_cohort2.xlsx") %>% mutate(Unique_Sample_ID = as.factor(Unique_Sample_ID))

smoking <- read.xlsx("iScience_input_data/smoking_info.xlsx", 1) %>% mutate(Unique_Sample_ID = as.factor(Unique_Sample_ID))
smoking <- smoking %>% mutate(Diagnosis_Smoking_classification = factor(Diagnosis_Smoking_classification, levels = c("Non", "Ex", "Current", "COPD")))


# Fig 2B ---------
## PCA lung
data <- elisa.lung.25.full %>% select(!starts_with("NEFA")) %>% 
  pivot_longer(where(is.numeric)) %>%
  mutate(Tissue = case_when(grepl("_S", name) ~ "Serum",
                            grepl("_L", name) ~ "Lung")) %>%
  mutate(name = gsub("_L", "", name),
         name = gsub("_S", "", name)) %>%
  pivot_wider()


PCA2(data, "Lung Combined cohort")

# select PC with the best diagnosis separation
axis.x = 2 # PC2 vs PC4 gives nice separation in the lung
axis.y = 4

p2l.combined<- fviz_pca_biplot(df1, title = NULL, 
                               axes = c(axis.x,axis.y), 
                               habillage= PC$Diagnosis, #label = "var",
                               palette = c("Donor" = "#757b87", "COPD" = "#791812"),
                               mean.point = FALSE, pointsize = 4, pointshape = 19, 
                               label ="var", select.var = list(contrib = 10), #addEllipses=TRUE, ellipse.level=0.95
                               labelsize = 4, arrowsize = 0.1, col.var = "#252525", repel = TRUE,
                               subtitle = paste0("Lung combined cohort", " PCA - Biplot")
)+
  expand_limits(x = c(min(df1$x[,paste0("PC", axis.x)])*1.1, max(df1$x[,axis.x])*1.1),
                y = c(min(df1$x[,paste0("PC", axis.y)])*1.1, max(df1$x[,axis.y])*1.1)) +  
  labs(x = paste0("PC", axis.x, " (", round(summary(df1)$importance[2,axis.x], digits = 3)*100, "%)"),
       y = paste0("PC", axis.y, " (", round(summary(df1)$importance[2,axis.y], digits = 3)*100, "%)"))+
  theme_pubr(base_size = 16)+
  
  theme(legend.position = "none") + theme_journal2


p2l.combined


# Fig 2C-E ------------
## Plasma
data <- elisa.serum.25.full %>% select(!starts_with("NEFA")) %>%
  pivot_longer(where(is.numeric)) %>%
  mutate(Tissue = case_when(grepl("_S", name) ~ "Serum",
                            grepl("_L", name) ~ "Lung")) %>%
  mutate(name = gsub("_L", "", name),
         name = gsub("_S", "", name)) %>%
  pivot_wider()


PCA2(data, "Plasma Combined cohort")

# select PC with the best diagnosis separation
axis.x = 1 # PC1 vs PC3 gives nice separation
axis.y = 3

p2p.combined<- fviz_pca_biplot(df1, title = NULL, 
                               axes = c(axis.x,axis.y), 
                               habillage= PC$Diagnosis, #label = "var",
                               palette = c("Donor" = "#757b87", "COPD" = "#791812"),
                               mean.point = FALSE, pointsize = 4, pointshape = 19, 
                               label ="var", select.var = list(contrib = 10), #addEllipses=TRUE, ellipse.level=0.95
                               labelsize = 4, arrowsize = 0.1, col.var = "#252525", repel = TRUE,
                               subtitle = paste0("Plasma combined cohort", " PCA - Biplot")
)+
  expand_limits(x = c(min(df1$x[,paste0("PC", axis.x)])*1.1, max(df1$x[,axis.x])*1.1),
                y = c(min(df1$x[,paste0("PC", axis.y)])*1.1, max(df1$x[,axis.y])*1.1)) +  
  labs(x = paste0("PC", axis.x, " (", round(summary(df1)$importance[2,axis.x], digits = 3)*100, "%)"),
       y = paste0("PC", axis.y, " (", round(summary(df1)$importance[2,axis.y], digits = 3)*100, "%)"))+
  theme_pubr(base_size = 16)+
    theme(legend.position = "none") + theme_journal2

p2p.combined



### PCA + smoking status
temp <- smoking %>% right_join(select(PC, -Diagnosis))

p4 <- temp %>%
  mutate(Smoking_py = as.numeric(Smoking_py)) %>% 
  ggplot(., aes(x= scale(PC1, center = F, scale = F), y=scale(PC3, center = F, scale = F), col=Smoking_py))+
  geom_point(data = subset(temp, !is.na(Smoking_py)), size=5, col="black")+
  geom_point(size=4)+
  labs(x = paste0("PC1 (", varPC1*100, "%)"), 
       y = paste0("PC3 (", varPC3*100, "%)"), 
       subtitle = "Smoking history")+
  theme(legend.position = "bottom")+
  
  theme(panel.grid.major = element_line(size = 0.75))+
  guides(col= guide_legend(title= "Pack years"))+
  
  expand_limits(x = c(min(df1$x[,axis.x])*1.1, max(df1$x[,axis.x])*1.1),
                y = c(min(df1$x[,axis.y])*1.1, max(df1$x[,axis.y])*1.1)) +  
  labs(x = paste0("PC", axis.x, " (", round(summary(df1)$importance[2,1], digits = 3)*100, "%)"),
       y = paste0("PC", axis.y, " (", round(summary(df1)$importance[2,3], digits = 3)*100, "%)"))+
  theme_journal2

p4


p3 <- temp %>%
  ggplot(., aes(x= scale(PC1, center = F, scale = F), y=scale(PC3, center = F, scale = F), col=Diagnosis_Smoking_classification))+
  
  geom_point(data = subset(temp, Diagnosis %in% c("Control")), size=5, col="black")+
  geom_point(size=4)+
  scale_color_manual(values = c("Non" = "#cccccc",
                                "Ex" = "#737373",
                                "Current" = "#252525",
                                "COPD" = "#791812"))+
  labs(x = paste0("PC1 (", varPC1*100, "%)"), 
       y = paste0("PC3 (", varPC3*100, "%)"), 
       subtitle = "Control smoking classification")+
  
  theme(legend.position = "bottom")+
  theme(panel.grid.major = element_line(size = 0.75))+
  guides(col= guide_legend(title= "Status"))+
  
  expand_limits(x = c(min(df1$x[,axis.x])*1.1, max(df1$x[,axis.x])*1.1),
                y = c(min(df1$x[,axis.y])*1.1, max(df1$x[,axis.y])*1.1)) +  
  labs(x = paste0("PC", axis.x, " (", round(summary(df1)$importance[2,1], digits = 3)*100, "%)"),
       y = paste0("PC", axis.y, " (", round(summary(df1)$importance[2,3], digits = 3)*100, "%)"))+
  theme_journal2

p3


### Patch figure ----------
patch <- p2l.combined + p2p.combined+p4+p3 +plot_layout(ncol = 4) #+plot_layout(axes = "collect")
patch

save_plot(paste0("plots/", Sys.Date(), "_ELISAS_combined_PCA_graphs_smoking", ".png"), plot = patch, 
          base_height = 6, base_width = 6.5*3.25, dpi = 600)


# Fig 2F --------
## Lung vs plasma biplot combined cohort
### sig.df
sig.df <- elisa.lung.25.full %>% 
  full_join(elisa.serum.25.full) %>%
  pivot_longer(where(is.numeric)) %>%
  mutate(Tissue = case_when(grepl("_S", name) ~ "Serum",
                            grepl("_L", name) ~ "Lung")) %>%
  mutate(name = gsub("_L", "", name),
         name = gsub("_S", "", name)) %>%
  mutate(value = log10(value)) %>%
  mutate(Diagnosis2 = gsub("Control", "Donor", Diagnosis)) %>%
  group_by(Tissue, name) %>% #
  summarise(pval = wilcox.test(value ~ Diagnosis2)$p.value) %>% 
  mutate(padj = p.adjust(pval, method = "BH", n = length(pval))) %>%
  mutate(signif = if_else(padj < 0.05, "ja", "nein")) %>%
  mutate(padj = -log10(padj)) %>%
  mutate(across(where(is.numeric), round, 5)) %>%
  select(-pval, -signif) %>%
  pivot_wider(names_from = Tissue, values_from = padj)

## Plasma lung plot
p <- elisa.lung.25.full %>% 
  full_join(elisa.serum.25.full) %>%
  pivot_longer(where(is.numeric)) %>%
  na.omit() %>%
  mutate(Tissue = case_when(grepl("_S", name) ~ "Serum",
                            grepl("_L", name) ~ "Lung")) %>%
  mutate(name = gsub("_L", "", name)) %>%
  mutate(name = gsub("_S", "", name)) %>%
  
  mutate(Diagnosis2 = gsub("Control", "Donor", Diagnosis)) %>%
  pivot_wider(names_from = name, values_from = value, values_fn = mean) %>% 
  group_by(Diagnosis2, Tissue) %>%
  summarise(across(where(is.numeric), mean, na.rm = T)) %>%
  
  pivot_longer(where(is.numeric))  %>%
  pivot_wider(names_from = Diagnosis2, values_from = value) %>% 
  
  mutate(Ave = rowMeans(select(., COPD, Donor))) %>%
  mutate(Diff = log2(COPD / Donor)) %>% #changed to L2FC
  
  select(Tissue, name, Diff) %>%
  pivot_wider(names_from = Tissue, values_from = Diff) %>%
  
  left_join(sig.df, by="name") %>% #sig.df contains info from combined cohort
  
  mutate(Significance = case_when(Lung.y <= 1.30103 & Serum.y <= 1.30103 ~ "p>0.05",
                                  Lung.y > 1.30103 & Serum.y <= 1.30103 ~ "Significant Lung",
                                  Lung.y <= 1.30103 & Serum.y > 1.30103 ~ "Significant Plasma",
                                  Lung.y > 1.30103 & Serum.y > 1.30103 ~ "Significant Lung + Plasma"
  )) %>% 
  
  
  na.omit() %>%
  
  rename(Parameter_Name = name) %>%
  left_join(select(Parameter_labels, Parameter_Name, expression)) %>%
  
  ggplot(., aes(x= Lung.x, y= Serum.x, col=Significance, label = expression))+
  geom_point(size=4,alpha=1)+
  
  geom_text_repel(col="black", parse = T) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+

  labs(x = expression(Lung~(Log["2"]~fold~change)), 
       y = expression(Plasma~(Log["2"]~fold~change)), 
       subtitle = "Compartmental differences in cytokine levels") +
  scale_x_continuous(limits = c(-2.75, 3.75))+
  scale_y_continuous(limits = c(-3, NA))+ 
  scale_color_manual(values = c("p>0.05" = "grey",
                                "Significant Lung" = "indianred3", 
                                "Significant Plasma" = "skyblue3",
                                "Significant Lung + Plasma" = "#673770"))+
  theme(panel.grid.major.y = element_blank())+
  theme(panel.grid.minor.y = element_blank())

plot(p)

plot(p)+ guides(col = guide_legend(nrow = 2))
p <- p + theme(legend.position = "none")

save_plot(paste0("plots/", Sys.Date(), "_", "ELISA_combined Serum vs Lung-L2FCv2_padj", ".png"), 
          plot = ggplot2::last_plot(), base_height = 7, base_width = 7, dpi = 600) 


# Fig 2G -----------
## Selected cytokines
factors <- c("IL6", "CCL5", "IL8", "CXCL9" ) 

elisa.25 <- elisa.serum.25.full %>% 
  full_join(elisa.lung.25.full) %>%
  pivot_longer(where(is.numeric)) %>%
  mutate(Tissue = case_when(grepl("_S", name) ~ "Serum",
                            grepl("_L", name) ~ "Lung")) %>%
  mutate(name = gsub("_L", "", name)) %>%
  mutate(name = gsub("_S", "", name)) %>%
  mutate(Diagnosis2 = gsub("Control", "Donor", Diagnosis), 
         Tissue.diagnosis = paste0(Tissue, ".", Diagnosis2)) %>%
  mutate(Tissue.diagnosis = factor(Tissue.diagnosis, levels = c("Lung.Donor", "Lung.COPD", "Serum.Donor", "Serum.COPD"))) %>%
  na.omit() %>%
  filter(name %in% factors) %>%
  mutate(value = log10(value)) %>%
  rename(Parameter_Name = name) %>%
  left_join(select(Parameter_labels, Parameter_Name, expression)) %>%
  mutate(Parameter_Name = factor(Parameter_Name, factors)) %>%

  ggplot(., aes(x = Tissue.diagnosis, y = value, fill = Diagnosis)) +
  geom_violin(color="white", alpha=0.4, lwd= 1)+
  geom_dotplot(dotsize = 1.5, binaxis="y", stackdir="center", stackgroups = TRUE,  binpositions = "all") +
  geom_vline(xintercept = 2.5, linetype = "dashed")+ 
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5)+
  labs(x = "", title = NULL, subtitle = NULL)+
  ylab(expression(paste("Concentration (LOG", " transformed)")))+
  stat_compare_means(comparisons = list( c("Serum.Donor", "Serum.COPD"), c("Lung.Donor", "Lung.COPD")), size = 6 ,
                     hide.ns = F, symnum.args = symnum.args,  
                     aes(label = paste0(after_stat(p.signif))))+ 
  scale_x_discrete(labels = c(Serum.Donor = "Donor \nPlasma", Serum.COPD = "COPD \nPlasma", Lung.Donor = "Donor \nLung", Lung.COPD = "COPD \nLung"))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  man.fill+
  axis+
  theme(legend.position="none")+
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
  facet_wrap(~Parameter_Name, scales = "free_y", 
             labeller=as_labeller(parameter_labeller, label_parsed), ncol = 2)

elisa.25

save_plot(paste0("plots/", Sys.Date(), "_ELISA_combined_dotplot_examples", ".png"),
          plot = elisa.25, 
          base_height = 8, base_width = 9, dpi = 300) 

# Fig 2I -----------
## see Fig2I_network.R script


# Fig S8B  ----------------
##plot remaining sig cytokines 
x <- sig.df %>% filter(Lung > -log10(0.05) | Serum > -log10(0.05)) %>% pull(name)

elisa.25 <- elisa.serum.25.full %>% 
  full_join(elisa.lung.25.full) %>%
  pivot_longer(where(is.numeric)) %>%
  mutate(Tissue = case_when(grepl("_S", name) ~ "Serum",
                            grepl("_L", name) ~ "Lung")) %>%
  mutate(name = gsub("_L", "", name)) %>%
  mutate(name = gsub("_S", "", name)) %>%
  mutate(Diagnosis2 = gsub("Control", "Donor", Diagnosis), 
         Tissue.diagnosis = paste0(Tissue, ".", Diagnosis2)) %>%
  mutate(Tissue.diagnosis = factor(Tissue.diagnosis, levels = c("Lung.Donor", "Lung.COPD", "Serum.Donor", "Serum.COPD"))) %>%
  na.omit() %>%
  
  filter(name %in% x) %>% 
  filter(!name %in% factors) %>%
  
  mutate(value = log10(value)) %>%
  rename(Parameter_Name = name) %>%
  left_join(select(Parameter_labels, Parameter_Name, expression)) %>%

  ggplot(., aes(x = Tissue.diagnosis, y = value, fill = Diagnosis)) +
  geom_violin(color="white", alpha=0.2, lwd= 1)+
  geom_dotplot(dotsize = 1.25, binaxis="y", stackdir="center", stackgroups = TRUE,  binpositions = "all") +
  
  geom_vline(xintercept = 2.5, linetype = "dashed")+ 
  
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5)+
  labs(x = "",
       title = "",
       subtitle = "Compartmental differences in cytokines levels")+
  ylab(expression(paste("Concentration (LOG", " transformed)")))+
  stat_compare_means(comparisons = list( c("Serum.Donor", "Serum.COPD"), c("Lung.Donor", "Lung.COPD")), size = 6 ,
                     hide.ns = F,  symnum.args = symnum.args,  aes(label = paste0(after_stat(p.signif))))+ 
  scale_x_discrete(labels = c(Serum.Donor = "Donor \nPlasma", Serum.COPD = "COPD \nPlasma", Lung.Donor = "Donor \nLung", Lung.COPD = "COPD \nLung"))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  man.fill+
  man.col+
  axis+
  theme(legend.position="none")+
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
  facet_wrap(~Parameter_Name, scales = "free_y", labeller=as_labeller(parameter_labeller, label_parsed), ncol = 4)

elisa.25

save_plot(paste0("plots/", Sys.Date(), "_ELISA_combined_remaining_dotplots", ".png"), plot = elisa.25, 
          base_height = 12, base_width = 16.5, dpi = 300) #ncol = 3 width = 13




# Fig S8A  -----------
sig.df <- 
  elisa.serum.25.full %>% left_join(smoking) %>% 
  pivot_longer(ends_with("_S")) %>%
  mutate(name = gsub("_S", "", name)) %>%
  mutate(value = log10(value)) %>%
  group_by(name) %>% 
  summarise(pval = kruskal.test(value ~ Diagnosis_Smoking_classification)$p.value) %>%
  mutate(padj = p.adjust(pval, method = "BH", n = length(pval))) %>%
  mutate(signif = case_when(padj <= 0.05 & padj > 0.01 ~"*", padj <= 0.01 & padj > 0.001 ~"**", padj <= 0.001 ~"***", .default = "ns"))

sig.list <- sig.df %>% filter(signif != "ns") %>% pull(name)


elisa.serum.25.full %>% left_join(smoking) %>% 
  pivot_longer(ends_with("_S")) %>%
  mutate(name = gsub("_S", "", name)) %>%
  mutate(value = log10(value)) %>%

  filter(name %in% sig.list) %>%
  
  ggplot(., aes(x = Diagnosis_Smoking_classification, y = value, fill = Diagnosis)) +
  geom_violin(color="white", alpha=0.2, lwd= 1)+
  geom_dotplot(dotsize = 2, binaxis="y", stackdir="center", stackgroups = TRUE,  binpositions = "all") +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5)+
  
    labs(x = "",
       title = "",
       subtitle = "Differences in plasma cytokines levels according to smoking status")+
  ylab(expression(paste("Concentration (LOG", " transformed)")))+
  stat_compare_means(comparisons = list( c("Non", "Ex"), c("Ex", "Current"), c("Non", "Current"), c("Non", "COPD"),
                                         c("Ex", "COPD"), c("Current", "COPD")),
                     hide.ns = TRUE, size = 6, symnum.args = symnum.args, aes(label = paste0(after_stat(p.signif))))+ 
  
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  man.fill+
  man.col+
  axis+
  theme(legend.position="none")+
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
  facet_wrap(~name, scales = "free_y",labeller=as_labeller(parameter_labeller, label_parsed), ncol = 5)

save_plot(paste0("plots/", Sys.Date(), "_ELISA_combined_sig_cytokines_smoking", ".png"), plot = ggplot2::last_plot(), 
          base_height = 12, base_width = 18, dpi = 600) #ncol = 3 width = 13



# Fig S6B -----------
## Lung FACS cohort 1
data <- elisa.lung.25.c1 %>% 
  pivot_longer(where(is.numeric)) %>%
  mutate(Tissue = case_when(grepl("_S", name) ~ "Serum",
                            grepl("_L", name) ~ "Lung")) %>%
  mutate(name = gsub("_L", "", name),
         name = gsub("_S", "", name)) %>%
  pivot_wider()

PCA2(data, "ELISA Lung FACS cohort 1")



# select PC with the best diagnosis separation
axis.x = 2
axis.y = 4
p2a <- fviz_pca_biplot(df1, title = NULL, 
                       axes = c(axis.x,axis.y), 
                       habillage= PC$Diagnosis, 
                       palette = c("Donor" = "#757b87", "COPD" = "#791812"),
                       mean.point = FALSE, pointsize = 4, pointshape = 19, 
                       label ="var", select.var = list(contrib = 10), 
                       labelsize = 4, arrowsize = 0.1, col.var = "#252525", repel = TRUE,
                       subtitle = paste0("Lung FACS cohort", " PCA - Biplot")
)+
  expand_limits(x = c(min(df1$x[,paste0("PC", axis.x)])*1.1, max(df1$x[,axis.x])*1.1),
                y = c(min(df1$x[,paste0("PC", axis.y)])*1.1, max(df1$x[,axis.y])*1.1)) +  
  labs(x = paste0("PC", axis.x, " (", round(summary(df1)$importance[2,axis.x], digits = 3)*100, "%)"),
       y = paste0("PC", axis.y, " (", round(summary(df1)$importance[2,axis.y], digits = 3)*100, "%)"))+
  theme_pubr(base_size = 16)+
  
  theme(legend.position = "none") + theme_journal2

p2a


# Fig S6C -----------
## Plasma FACS cohort 1 
data <- elisa.serum.25.c1 %>% #filter(Unique_Sample_ID %in% list.facs.s) %>% select(!starts_with("NEFA")) %>%
  pivot_longer(where(is.numeric)) %>%
  mutate(Tissue = case_when(grepl("_S", name) ~ "Serum",
                            grepl("_L", name) ~ "Lung")) %>%
  mutate(name = gsub("_L", "", name),
         name = gsub("_S", "", name)) %>%
  pivot_wider()

PCA2(data, "Plasma FACS cohort")

# select PC with the best diagnosis separation
axis.x = 2
axis.y = 3
p2b <- fviz_pca_biplot(df1, title = NULL, 
                       axes = c(axis.x,axis.y), 
                       habillage= PC$Diagnosis, #label = "var",
                       palette = c("#757b87", "#791812"),
                       mean.point = FALSE, pointsize = 4, pointshape = 19, 
                       label ="var", select.var = list(contrib = 10), #addEllipses=TRUE, ellipse.level=0.95
                       labelsize = 4, arrowsize = 0.1, col.var = "#252525", repel = TRUE,
                       subtitle = paste0("Plasma FACS cohort", " PCA - Biplot")
)+
  expand_limits(x = c(min(df1$x[,paste0("PC", axis.x)])*1.1, max(df1$x[,axis.x])*1.1),
                y = c(min(df1$x[,paste0("PC", axis.y)])*1.1, max(df1$x[,axis.y])*1.1)) +  
  labs(x = paste0("PC", axis.x, " (", round(summary(df1)$importance[2,axis.x], digits = 3)*100, "%)"),
       y = paste0("PC", axis.y, " (", round(summary(df1)$importance[2,axis.y], digits = 3)*100, "%)"))+
  theme_pubr(base_size = 16)+
  
  theme(legend.position = "none") + theme_journal2

p2b



# Fig S6D ---------------------------
### sig.df 
sig.df <- elisa.lung.25.c1 %>% 
  full_join(elisa.serum.25.c1) %>%
  pivot_longer(where(is.numeric)) %>%
  mutate(Tissue = case_when(grepl("_S", name) ~ "Serum",
                            grepl("_L", name) ~ "Lung")) %>%
  mutate(name = gsub("_L", "", name),
         name = gsub("_S", "", name)) %>%
  mutate(value = log10(value)) %>%
  
  mutate(Diagnosis2 = gsub("Control", "Donor", Diagnosis)) %>%
  group_by(Tissue, name) %>% 
  summarise(pval = wilcox.test(value ~ Diagnosis2)$p.value) %>% 
  mutate(padj = p.adjust(pval, method = "BH", n = length(pval))) %>%
  mutate(padj = -log10(padj)) %>%
  mutate(signif = if_else(pval < 0.05, "ja", "nein")) %>%
  select(-pval, -signif) %>%
  pivot_wider(names_from = Tissue, values_from = padj)


## Plasma lung plot 
p1 <- elisa.lung.25.c1 %>% 
  full_join(elisa.serum.25.c1) %>%
  
  pivot_longer(where(is.numeric)) %>%
  na.omit() %>%
  mutate(Tissue = case_when(grepl("_S", name) ~ "Serum",
                            grepl("_L", name) ~ "Lung")) %>%
  mutate(name = gsub("_L", "", name)) %>%
  mutate(name = gsub("_S", "", name)) %>%
  
  #filter(!name %in% c("NEFA", "CCL5", "TSLP")) %>% #NEFA removed !
  
  mutate(Diagnosis2 = gsub("Control", "Donor", .$Diagnosis)) %>%
  pivot_wider(names_from = name, values_from = value, values_fn = mean) %>% 
  group_by(Diagnosis2, Tissue) %>%
  summarise(across(where(is.numeric), mean, na.rm = T)) %>%
  
  pivot_longer(where(is.numeric))  %>%
  pivot_wider(names_from = Diagnosis2, values_from = value) %>% 
  mutate(Ave = rowMeans(select(., COPD, Donor))) %>%
  mutate(Diff = log2(COPD / Donor)) %>% 

  select(Tissue, name, Diff) %>%
  pivot_wider(names_from = Tissue, values_from = Diff) %>%
  
  left_join(sig.df, by="name") %>% 
  
  mutate(Significance = case_when(Lung.y <= 1.30103 & Serum.y <= 1.30103 ~ "p>0.05",
                                  Lung.y > 1.30103 & Serum.y <= 1.30103 ~ "Significant Lung",
                                  Lung.y <= 1.30103 & Serum.y > 1.30103 ~ "Significant Plasma",
                                  Lung.y > 1.30103 & Serum.y > 1.30103 ~ "Significant Lung + Plasma"
  )) %>% 
  
  
  na.omit() %>%
  
  rename(Parameter_Name = name) %>%
  left_join(select(Parameter_labels, Parameter_Name, expression)) %>%
  
  ggplot(., aes(x= Lung.x, y= Serum.x, col=Significance, label = expression))+
  geom_point(size=4,alpha=1)+
  geom_text_repel(col="black", parse = T) +

  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  

  labs(x = expression(Lung~(Log["2"]~fold~change)), 
       y = expression(Plasma~(Log["2"]~fold~change)), 
       subtitle = "ELISA cohort1") +
  scale_x_continuous(limits = c(-2.75, 3.75))+ 
  scale_y_continuous(limits = c(-3, NA))+ 
  scale_color_manual(values = c("p>0.05" = "grey",
                                "Significant Lung" = "indianred3", 
                                "Significant Plasma" = "skyblue3",
                                "Significant Lung + Plasma" = "#673770"))+
  theme(panel.grid.major.y = element_blank())+
  theme(panel.grid.minor.y = element_blank())

plot(p1)

plot(p1)+ guides(col = guide_legend(nrow = 2))
p1.facs <- p1 + theme(legend.position = "none")+labs(subtitle = "")



## Patch figure
patch <- p2a+theme(legend.position = "none")+ p2b+theme(legend.position = "none")+
  p1.facs
patch

save_plot(paste0("plots/", Sys.Date(), "_FACS_cohort1 cytokine overview", ".png"), plot = patch, 
          base_height = 6, base_width = 16, dpi = 600)



# Fig S6E ----------------
## Selected cytokines
factors <- c("IL6", "GM_CSF", "IFN_la1", "CXCL9") 

elisa.25 <- elisa.serum.25.c1 %>% 
  full_join(elisa.lung.25.c1) %>%
  pivot_longer(where(is.numeric)) %>%
  mutate(Tissue = case_when(grepl("_S", name) ~ "Serum",
                            grepl("_L", name) ~ "Lung")) %>%
  mutate(name = gsub("_L", "", name)) %>%
  mutate(name = gsub("_S", "", name)) %>%
  mutate(Diagnosis2 = gsub("Control", "Donor", Diagnosis), 
         Tissue.diagnosis = paste0(Tissue, ".", Diagnosis2)) %>%
  mutate(Tissue.diagnosis = factor(Tissue.diagnosis, levels = c("Lung.Donor", "Lung.COPD", "Serum.Donor", "Serum.COPD"))) %>%
  na.omit() %>%
  filter(name %in% factors) %>%
  mutate(value = log10(value)) %>%
  rename(Parameter_Name = name) %>%
  left_join(select(Parameter_labels, Parameter_Name, expression)) %>%
  mutate(Parameter_Name = factor(Parameter_Name, factors)) %>%
  
  ggplot(., aes(x = Tissue.diagnosis, y = value, fill = Diagnosis)) +
  geom_violin(color="white", alpha=0.4, lwd= 1)+
  geom_dotplot(dotsize = 1.5, binaxis="y", stackdir="center", stackgroups = TRUE,  binpositions = "all") +
  geom_vline(xintercept = 2.5, linetype = "dashed")+ 
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5)+
  labs(x = "", title = NULL, subtitle = NULL)+
  ylab(expression(paste("Concentration (LOG", " transformed)")))+
  stat_compare_means(comparisons = list( c("Serum.Donor", "Serum.COPD"), c("Lung.Donor", "Lung.COPD")), size = 6 ,
                     hide.ns = F, symnum.args = symnum.args,  
                     aes(label = paste0(after_stat(p.signif))))+ 
  scale_x_discrete(labels = c(Serum.Donor = "Donor \nPlasma", Serum.COPD = "COPD \nPlasma", Lung.Donor = "Donor \nLung", Lung.COPD = "COPD \nLung"))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  man.fill+
  axis+
  theme(legend.position="none")+
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
  facet_wrap(~Parameter_Name, scales = "free_y", 
             labeller=as_labeller(parameter_labeller, label_parsed), ncol = 4)
elisa.25

save_plot(paste0("plots/", Sys.Date(), "_ELISA_cohort1_dotplot_examples", ".png"), plot = elisa.25, 
          base_height = 5, base_width = 16, dpi = 300) #ncol = 3 width = 13, base_height = 8, base_width = 8.5


# Figure S6 extra ----------------------------
### FACS cohort 1 all dotblots expect those in main figure
elisa.lung.25.c1 %>% 
  full_join(elisa.serum.25.c1) %>% mutate(FACS = "FACS") %>%
  pivot_longer(where(is.numeric)) %>%
  mutate(Tissue = case_when(grepl("_S", name) ~ "Serum",
                            grepl("_L", name) ~ "Lung")) %>%
  mutate(name = gsub("_L", "", name)) %>%
  mutate(name = gsub("_S", "", name)) %>%
  mutate(Diagnosis2 = gsub("Control", "Donor", Diagnosis), 
         Tissue.diagnosis = paste0(Tissue, ".", Diagnosis2)) %>%
  mutate(Tissue.diagnosis = factor(Tissue.diagnosis, levels = c("Lung.Donor", "Lung.COPD", "Serum.Donor", "Serum.COPD"))) %>%
  na.omit() %>%
  mutate(value = log10(value)) %>%
  filter(!name %in% factors) %>%
  
  rename(Parameter_Name = name) %>%
  left_join(select(Parameter_labels, Parameter_Name, expression)) %>%
  ggplot(., aes(x = Tissue.diagnosis, y = value, fill = Diagnosis)) +
  geom_violin(color="white", alpha=0.2, lwd= 1)+
  geom_dotplot(dotsize = .75, binaxis="y", stackdir="center", stackgroups = TRUE,  binpositions = "all") +
  
  geom_vline(xintercept = 2.5, linetype = "dashed")+ 
  
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5)+
  labs(x = "", y = expression(paste("Concentration (Log" [10], " transformed)")),
       subtitle = "Compartmental differences in cytokines levels")+
  
    stat_compare_means(comparisons = list( c("Serum.Donor", "Serum.COPD"), c("Lung.Donor", "Lung.COPD")), hide.ns = T)+ 
  scale_x_discrete(labels = c(Serum.Donor = "Donor \nPlasma", Serum.COPD = "COPD \nPlasma", Lung.Donor = "Donor \nLung", Lung.COPD = "COPD \nLung"))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  man.fill+
  man.col+
  axis+
  theme(legend.position="none")+
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
  facet_wrap(~Parameter_Name, scales = "free_y", labeller=as_labeller(parameter_labeller, label_parsed))

save_plot(paste0("plots/", Sys.Date(), "_", "ELISA_cohort1_dotplot_all remaing cytokines", ".png"), 
          plot = ggplot2::last_plot(), base_height = 20, base_width = 27, dpi = 600)


# Fig S7B -----------
## Lung cohort 2
data <- elisa.lung.25.c2 %>% #filter(!Unique_Sample_ID %in% list.facs.l) %>% select(!starts_with("NEFA")) %>%
  pivot_longer(where(is.numeric)) %>%
  mutate(Tissue = case_when(grepl("_S", name) ~ "Serum",
                            grepl("_L", name) ~ "Lung")) %>%
  mutate(name = gsub("_L", "", name),
         name = gsub("_S", "", name)) %>%
  pivot_wider()

PCA2(data, "Lung, cohort 2")


# select PC with the best diagnosis separation
axis.x = 2
axis.y = 3
p2c <- fviz_pca_biplot(df1, title = NULL, 
                       axes = c(axis.x,axis.y), 
                       habillage= PC$Diagnosis, #label = "var",
                       palette = c("Donor" = "#757b87", "COPD" = "#791812"),
                       mean.point = FALSE, pointsize = 4, pointshape = 19, 
                       label ="var", select.var = list(contrib = 10), #addEllipses=TRUE, ellipse.level=0.95
                       labelsize = 4, arrowsize = 0.1, col.var = "#252525", repel = TRUE,
                       subtitle = paste0("Lung cohort 2", " PCA - Biplot")
)+
  expand_limits(x = c(min(df1$x[,paste0("PC", axis.x)])*1.1, max(df1$x[,axis.x])*1.1),
                y = c(min(df1$x[,paste0("PC", axis.y)])*1.1, max(df1$x[,axis.y])*1.1)) +  
  labs(x = paste0("PC", axis.x, " (", round(summary(df1)$importance[2,axis.x], digits = 3)*100, "%)"),
       y = paste0("PC", axis.y, " (", round(summary(df1)$importance[2,axis.y], digits = 3)*100, "%)"))+
  theme_pubr(base_size = 16)+
  
  theme(legend.position = "none") + theme_journal2

p2c


# Fig S7C -----------
## Plasma nonFACS cohort 
data <- elisa.serum.25.c2 %>%
  
  pivot_longer(where(is.numeric)) %>%
  mutate(Tissue = case_when(grepl("_S", name) ~ "Serum",
                            grepl("_L", name) ~ "Lung")) %>%
  mutate(name = gsub("_L", "", name),
         name = gsub("_S", "", name)) %>%
  pivot_wider()

PCA2(data, "Plasma, cohort 2")

p <- data %>% column_to_rownames(var = "Unique_Sample_ID") %>% select(where(is.numeric))
setdiff(names(p), names(p[y,x]))

# select PC with the best diagnosis separation
axis.x = 1
axis.y = 2
p2d <- fviz_pca_biplot(df1, title = NULL, 
                       axes = c(axis.x,axis.y), 
                       habillage= PC$Diagnosis, #label = "var",
                       palette = c("Donor" = "#757b87", "COPD" = "#791812"),
                       mean.point = FALSE, pointsize = 4, pointshape = 19, 
                       label ="var", select.var = list(contrib = 10), #addEllipses=TRUE, ellipse.level=0.95
                       labelsize = 4, arrowsize = 0.1, col.var = "#252525", repel = TRUE,
                       subtitle = paste0("Plasma cohort 2", " PCA - Biplot")
)+
  expand_limits(x = c(min(df1$x[,paste0("PC", axis.x)])*1.1, max(df1$x[,axis.x])*1.1),
                y = c(min(df1$x[,paste0("PC", axis.y)])*1.1, max(df1$x[,axis.y])*1.1)) +  
  labs(x = paste0("PC", axis.x, " (", round(summary(df1)$importance[2,axis.x], digits = 3)*100, "%)"),
       y = paste0("PC", axis.y, " (", round(summary(df1)$importance[2,axis.y], digits = 3)*100, "%)"))+
  theme_pubr(base_size = 16)+
  
  theme(legend.position = "none") + theme_journal2


p2d

# Fig S7D ----------------
## Plasma lung plot
## sig.df
sig.df <- elisa.lung.25.c2 %>% 
    full_join(elisa.serum.25.c2) %>%
  pivot_longer(where(is.numeric)) %>%
  mutate(Tissue = case_when(grepl("_S", name) ~ "Serum",
                            grepl("_L", name) ~ "Lung")) %>%
  mutate(name = gsub("_L", "", name),
         name = gsub("_S", "", name)) %>%
  mutate(value = log10(value)) %>%
  
  mutate(Diagnosis2 = gsub("Control", "Donor", Diagnosis)) %>%
  group_by(Tissue, name) %>% 
  summarise(pval = wilcox.test(value ~ Diagnosis2)$p.value) %>% 
  mutate(padj = p.adjust(pval, method = "BH", n = length(pval))) %>%
  mutate(padj = -log10(padj)) %>%
  mutate(signif = if_else(pval < 0.05, "ja", "nein")) %>%
  select(-pval, -signif) %>%
  pivot_wider(names_from = Tissue, values_from = padj)



## Plasma lung plot
p1 <- elisa.lung.25.c2 %>% 
  full_join(elisa.serum.25.c2) %>%
  
  pivot_longer(where(is.numeric)) %>%
  na.omit() %>%
  mutate(Tissue = case_when(grepl("_S", name) ~ "Serum",
                            grepl("_L", name) ~ "Lung")) %>%
  mutate(name = gsub("_L", "", name)) %>%
  mutate(name = gsub("_S", "", name)) %>%
  
  mutate(Diagnosis2 = gsub("Control", "Donor", .$Diagnosis)) %>%
  pivot_wider(names_from = name, values_from = value, values_fn = mean) %>% 
  group_by(Diagnosis2, Tissue) %>%
  summarise(across(where(is.numeric), mean, na.rm = T)) %>%
  
  pivot_longer(where(is.numeric))  %>%
  pivot_wider(names_from = Diagnosis2, values_from = value) %>% 
  
  mutate(Ave = rowMeans(select(., COPD, Donor))) %>%
  mutate(Diff = log2(COPD / Donor)) %>% 
  
  select(Tissue, name, Diff) %>%
  pivot_wider(names_from = Tissue, values_from = Diff) %>%
  
  left_join(sig.df, by="name") %>% 
  
  mutate(Significance = case_when(Lung.y <= 1.30103 & Serum.y <= 1.30103 ~ "p>0.05",
                                  Lung.y > 1.30103 & Serum.y <= 1.30103 ~ "Significant Lung",
                                  Lung.y <= 1.30103 & Serum.y > 1.30103 ~ "Significant Plasma",
                                  Lung.y > 1.30103 & Serum.y > 1.30103 ~ "Significant Lung + Plasma"
  )) %>% 
  
  
  na.omit() %>%
  
  rename(Parameter_Name = name) %>%
  left_join(select(Parameter_labels, Parameter_Name, expression)) %>%
  
  ggplot(., aes(x= Lung.x, y= Serum.x, col=Significance, label = expression))+
  geom_point(size=4,alpha=1)+
  
  geom_text_repel(col="black", parse = T) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  labs(x = expression(Lung~(Log["2"]~fold~change)), 
       y = expression(Plasma~(Log["2"]~fold~change)), 
       subtitle = "nonFACS cohort 2") +
  scale_x_continuous(limits = c(-2.75, 3.75))+ 
  scale_y_continuous(limits = c(-3, NA))+ 
  scale_color_manual(values = c("p>0.05" = "grey",
                                "Significant Lung" = "indianred3", 
                                "Significant Plasma" = "skyblue3",
                                "Significant Lung + Plasma" = "#673770"))+
  theme(panel.grid.major.y = element_blank())+
  theme(panel.grid.minor.y = element_blank())

plot(p1)

plot(p1)+ guides(col = guide_legend(nrow = 2))
p1.nonfacs <- p1 + theme(legend.position = "none")+ labs(subtitle = "")



### Patch figure ------
patch <- p2c+theme(legend.position = "none")+ p2d+theme(legend.position = "none")+
  p1.nonfacs+theme(legend.position = "none")

patch

save_plot(paste0("plots/", Sys.Date(), "_cohort 2 cytokine overview", ".png"), plot = patch, 
          base_height = 6, base_width = 16, dpi = 600)


### Selected cytokines  ----------------
factors <- c("IL6", "CCL5", "CCL2", "CXCL9") 

elisa.25.non <- elisa.lung.25.c2 %>% 
  full_join(elisa.serum.25.c2) %>%
  pivot_longer(where(is.numeric)) %>%
  mutate(Tissue = case_when(grepl("_S", name) ~ "Serum",
                            grepl("_L", name) ~ "Lung")) %>%
  mutate(name = gsub("_L", "", name)) %>%
  mutate(name = gsub("_S", "", name)) %>%
  mutate(Diagnosis2 = gsub("Control", "Donor", Diagnosis), 
         Tissue.diagnosis = paste0(Tissue, ".", Diagnosis2)) %>%
  mutate(Tissue.diagnosis = factor(Tissue.diagnosis, levels = c("Lung.Donor", "Lung.COPD", "Serum.Donor", "Serum.COPD"))) %>%
  na.omit() %>%
  filter(name %in% factors) %>%
  mutate(value = log10(value)) %>%
  
  rename(Parameter_Name = name) %>%
  left_join(select(Parameter_labels, Parameter_Name, expression)) %>%
  mutate(Parameter_Name = factor(Parameter_Name, factors)) %>%
  
  
  ggplot(., aes(x = Tissue.diagnosis, y = value, fill = Diagnosis)) +
  geom_violin(color="white", alpha=0.2, lwd= 1)+
  geom_dotplot(dotsize = 1.5, binaxis="y", stackdir="center", stackgroups = TRUE,  binpositions = "all") +
  
  geom_vline(xintercept = 2.5, linetype = "dashed")+ 
  
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5)+
  labs(x = "",
       title = "",
       subtitle = "")+
  ylab(expression(paste("Concentration (LOG", " transformed)")))+
  stat_compare_means(comparisons = list( c("Serum.Donor", "Serum.COPD"), c("Lung.Donor", "Lung.COPD")), size = 6 ,
                     hide.ns = F, symnum.args = symnum.args,  
                     aes(label = paste0(after_stat(p.signif))))+ 
  scale_x_discrete(labels = c(Serum.Donor = "Donor \nPlasma", Serum.COPD = "COPD \nPlasma", Lung.Donor = "Donor \nLung", Lung.COPD = "COPD \nLung"))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  man.fill+
  man.col+
  axis+
  theme(legend.position="none")+
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
  facet_wrap(~Parameter_Name, scales = "free_y", labeller=as_labeller(parameter_labeller, label_parsed), ncol = 4)

elisa.25.non

save_plot(paste0("plots/", Sys.Date(), "_ELISA_cohort2_dotplot_examples", ".png"), 
          plot = elisa.25.non, 
          base_height = 5, base_width = 16, dpi = 300)


# Fig S7 extra -----------
## non-facs dotplots
elisa.lung.25.c2 %>% 
  full_join(elisa.serum.25.c2) %>% mutate(FACS = "non") %>%
  
  pivot_longer(where(is.numeric)) %>%
  mutate(Tissue = case_when(grepl("_S", name) ~ "Serum",
                            grepl("_L", name) ~ "Lung")) %>%
  mutate(name = gsub("_L", "", name)) %>%
  mutate(name = gsub("_S", "", name)) %>%
  mutate(Diagnosis2 = gsub("Control", "Donor", Diagnosis), 
         Tissue.diagnosis = paste0(Tissue, ".", Diagnosis2)) %>%
  mutate(Tissue.diagnosis = factor(Tissue.diagnosis, levels = c("Lung.Donor", "Lung.COPD", "Serum.Donor", "Serum.COPD"))) %>%
  na.omit() %>%
  filter(!name %in% factors) %>%
  mutate(value = log10(value)) %>%

  rename(Parameter_Name = name) %>%
  left_join(select(Parameter_labels, Parameter_Name, expression)) %>%

  ggplot(., aes(x = Tissue.diagnosis, y = value, fill = Diagnosis)) +
  geom_violin(color="white", alpha=0.2, lwd= 1)+
  geom_dotplot(dotsize = .75, binaxis="y", stackdir="center", stackgroups = TRUE,  binpositions = "all") +
  
  geom_vline(xintercept = 2.5, linetype = "dashed")+ 
  
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5)+
  labs(x = "",
       title = "cohort 2",
       subtitle = "Compartmental differences in cytokines levels")+
  ylab(expression(paste("Concentration (Log" [10], " transformed)")))+
  stat_compare_means(comparisons = list( c("Serum.Donor", "Serum.COPD"), c("Lung.Donor", "Lung.COPD")), hide.ns = T)+ 
  scale_x_discrete(labels = c(Serum.Donor = "Donor \nPlasma", Serum.COPD = "COPD \nPlasma", Lung.Donor = "Donor \nLung", Lung.COPD = "COPD \nLung"))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  man.fill+
  man.col+
  axis+
  theme(legend.position="none")+
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
  facet_wrap(~Parameter_Name, scales = "free_y", labeller=as_labeller(parameter_labeller, label_parsed))

save_plot(paste0("plots/", Sys.Date(), "_", "E123 Serum vs Lung-DotPlots_nonFACScohort", ".png"), 
          plot = ggplot2::last_plot(), base_height = 20, base_width = 27, dpi = 600)


