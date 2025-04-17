# Fig 4 ----------------

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
library(rstatix) #cohens n
library(readxl)
library(tidyr)
library(missMDA)
library(pheatmap)


## Formating ---------------
theme_set(theme_pubr(base_size=16, border = T))

cons_man.fill <- scale_fill_manual(values=c("A" = "#661510", "B" = "#e34e45")) #A darker B lighter
cons_man.col <- scale_color_manual(values=c("A" = "#661510", "B" = "#e34e45")) #A darker B lighter

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

# Import data -------------------
Parameter_labels <- read_excel("iScience_input_data/Parameter_labels.xlsx", 1)
parameter_labeller <- Parameter_labels$expression
names(parameter_labeller) <- Parameter_labels$Parameter_Name

facs.per.trans <- read.xlsx("iScience_input_data/hFACS_Donor_COPD_per_LOG_x1c_full.xlsx") %>% mutate(Unique_Sample_ID = as.factor(Unique_Sample_ID))
elisa.serum.25.c1 <- read.xlsx("iScience_input_data/hELISAS_Donor_COPD_conc_LOG_0quarterofmin_cohort1.xlsx") %>% mutate(Unique_Sample_ID = as.factor(Unique_Sample_ID))

COPD_subclass <- read.xlsx("iScience_input_data/COPD_subclass.xlsx") %>% 
  mutate(Unique_Sample_ID = as.factor(Lung.Matchcode)) %>%
  mutate(COPD_subclass = as.factor(COPD_subclass))

gas.exchange <- read.xlsx("iScience_input_data/gas_exchange.xlsx") %>% 
  mutate(Unique_Sample_ID = as.factor(Unique_Sample_ID)) %>%
  mutate(COPD_subclass = as.factor(COPD_subclass))



# Fig 4B ---------
## see following script ------------
#R_Fig4Btop_PCA_biplot.r
#R_Fig4_COPD_subclass.r
#R_Fig4Bdown_OPLSDA_Scores_plot.r

# Fig 4C -----------
## see following script ------------
#R_Fig4Ctop_PCA_biplot.r
#R_Fig4Cdown_OPLSDA_Scores_plot.r

# Fig 4D -----------
## see following script ------------
#R_Fig4D_Heatmap.r

# Fig 4E ---------------
## Cohen plot FACS+plasma
data <- COPD_subclass %>%
  left_join(facs.per.trans) %>% 
  select(!c("Lymphocytes", 
            "Macrophages", 
            "Mono_clas", "Mono_int", "Mono_non", 
            "PMNL",
            "Mast", 
            "Macs_CD14med", "Macs_CD14hi", # in first attempt Macs_CD14med was twice and no Macs_CD14hi
            "DC_CD209neg_CD11cpos", "DC_CD209pos_CD11cneg", "DC_CD209pos_CD11cpos", 
            "GRAN_CD193pos_CD16neg", "GRAN_CD193pos_CD16pos")) %>% 
  filter(Diagnosis != "Donor") %>%
  left_join((elisa.serum.25.c1 %>% mutate(across(where(is.numeric), log10)))) %>% 
  column_to_rownames(var = "Unique_Sample_ID")


### cohens_d
temp <- data %>% pivot_longer(where(is.numeric)) %>%
  group_by(name) %>% 
  cohens_d(value ~ COPD_subclass, paired = F) %>%
  select(name, effsize, magnitude)


### sig.df
sig.df <- data %>%
  pivot_longer(where(is.numeric)) %>%
  na.omit() %>%
  group_by(name) %>% 
  summarise(pval = wilcox.test(value ~ COPD_subclass)$p.value) %>% 
  mutate(padj = p.adjust(pval, method = "BH", n = length(pval))) %>%
  mutate(pval = round(pval, 5), padj = round(padj, 5)) %>%
  mutate(log.padj = round(-log10(padj), 5)) %>%
  mutate(signif = case_when(padj <= 0.05 & padj > 0.01 ~"*", padj <= 0.01 & padj > 0.001 ~"**", padj <= 0.001 ~"***", .default = "ns")) %>%
  left_join(temp)

sig.list <- sig.df %>% filter(pval < 0.05) %>% pull(name)

### graph 
sig.df %>% 
  mutate(name = gsub("_S", "", name)) %>%
  rename(Parameter_Name = name) %>%
  left_join(select(Parameter_labels, Parameter_Name, expression)) %>%
  arrange(desc(magnitude)) %>% 
  
  mutate(expression = fct_reorder(expression, abs(effsize))) %>% 
  na.omit() %>%
  
  ggplot(aes(effsize, expression, colour = magnitude))+
  geom_point(aes(size = -log10(padj)))+
  labs(x = "cohens_d", y="", title = "Effect size") +
  scale_x_continuous(limits = c(-1.5, NA))+
  scale_y_discrete(labels = scales::label_parse())+
  scale_color_manual(values = c("negligible" = "gray50", "small" = "#fee090", "moderate" = "#fc8d59", "large" = "#d73027"))+
  theme(legend.position = "right")+
  guides(colour = guide_legend(ncol= 1)) + theme_journal2


save_plot(paste0("plots/", Sys.Date(), "_COPD_subclass_FACS_Serum-effect-size", ".png"), 
          plot = ggplot2::last_plot(), base_height = 12, base_width = 10, dpi = 300)



# Fig 4F -----------------
## plot different cell types
data <- COPD_subclass %>%
  left_join(facs.per.trans) %>% 
  select(!c("Lymphocytes", 
            "Macrophages", 
            "Mono_clas", "Mono_int", "Mono_non", 
            "PMNL",
            "Mast", 
            "Macs_CD14med", "Macs_CD14hi", # in first attempt Macs_CD14med was twice and no Macs_CD14hi
            "DC_CD209neg_CD11cpos", "DC_CD209pos_CD11cneg", "DC_CD209pos_CD11cpos", 
            "GRAN_CD193pos_CD16neg", "GRAN_CD193pos_CD16pos")) %>% 
  filter(Diagnosis != "Donor") %>%
  column_to_rownames(var = "Unique_Sample_ID")

data %>%
  rownames_to_column(var = "Unique_Sample_ID") %>%
  pivot_longer(where(is.numeric)) %>%
  filter(name %in% sig.list) %>%
  
  mutate(name = factor(name, levels = c("CD8", "DC_CD209pos_CD11cposCD1a", 
                                        "Macs_CD14hi_CD1aposHLApos","CD45pos_CD203pos_CD117pos",
                                        "DC_CD209neg_CD11cposCd1a", "Macs_CD14med_CD1aposHLApos"))) %>%
  
  ggplot(., aes(x = COPD_subclass , y = value, fill = COPD_subclass ))+
  geom_violin(color="white", alpha=0.2, lwd= 1)+
  geom_dotplot(dotsize = 2, binaxis="y", stackdir="center") +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5)+
  labs(x = "",  y= expression("%"*CD45^"+"~cells~(LOG~transformed)))+
  stat_compare_means(aes(label = paste0(after_stat(p.signif))), size = 8,
                     symnum.args = symnum.args )+ 
  scale_fill_manual(values=c("A" = "#661510", "B" = "#e34e45", "Donor"  = "gray75"))+
  scale_x_discrete(labels=c("A" = "Subgroup I", "B" = "Subgroup II"))+
  axis+
  theme(legend.position="none")+
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.25)))+
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
  facet_wrap(~name, scales = "fixed", labeller=as_labeller(parameter_labeller, label_parsed), ncol = 3)
  


save_plot(paste0("plots/", Sys.Date(), "_COPD_subclass_FACS_Serum-FACS", ".png"), 
          plot = ggplot2::last_plot(), base_height = 10, base_width = 9.5, dpi = 600)

# Fig 4G -----------------
## plot different cytokines
data <- COPD_subclass %>%
  left_join(elisa.serum.25.c1) %>% 
  mutate(across(where(is.numeric), log10)) %>% 
  filter(Diagnosis != "Donor") %>%
  column_to_rownames(var = "Unique_Sample_ID")


data %>%
  rownames_to_column(var = "Unique_Sample_ID") %>%
  pivot_longer(where(is.numeric)) %>%
  filter(name %in% sig.list) %>%
  
  ggplot(., aes(x = COPD_subclass , y = value, fill = COPD_subclass ))+
  geom_violin(color="white", alpha=0.2, lwd= 1)+
  geom_dotplot(dotsize = 2, binaxis="y", stackdir="center") +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5)+
  labs(x = "", y = expression(paste("Concentration (LOG", " transformed)")))+
  stat_compare_means(aes(label = paste0(after_stat(p.signif))), size = 8,
                     symnum.args = symnum.args )+ 
  scale_fill_manual(values=c("A" = "#661510", "B" = "#e34e45", "Donor"  = "gray75"))+
  scale_x_discrete(labels=c("A" = "Subgroup I", "B" = "Subgroup II"))+
  axis+
  theme(legend.position="none")+
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.25)))+
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
  facet_wrap(~name, scales = "fixed", labeller=as_labeller(parameter_labeller, label_parsed), ncol = 5)

save_plot(paste0("plots/", Sys.Date(), "_COPD_subclass_FACS_Serum-Cytokines", ".png"), plot = ggplot2::last_plot(), 
          base_height = 5.5, base_width = 13, dpi = 300)


# Fig 5B --------------
f5b1 <- gas.exchange %>% 
  pivot_longer(where(is.numeric)) %>%
  filter(name %in% c("Airspace.enlargement")) %>% 
  mutate(name = gsub("Airspace.enlargement", "Airspace enlargement", name)) %>%
  
  ggplot(aes(COPD_subclass, value, fill = COPD_subclass))+
  geom_dotplot(dotsize = 1.3, binaxis="y", stackdir="center") +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5)+
  
  labs(x = "", y = "Mean interseptal distance \n(µm)")+
  stat_compare_means(aes(label = paste0(after_stat(p.signif))), size = 8,
                     symnum.args = symnum.args )+ 
  
  axis+
  theme(legend.position="none")+
  
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  scale_fill_manual(values=c("A" = "#661510", "B" = "#e34e45", "Donor"  = "gray75"))+
  scale_x_discrete(labels=c("A" = "Subgroup I", "B" = "Subgroup II"))+
  
  facet_wrap(~name, scales = "free_y", ncol = 2)

f5b1


f5b2 <- gas.exchange %>% 
  pivot_longer(where(is.numeric)) %>%
  filter(name %in% c("pO2_mmHg")) %>% 
  mutate(name = gsub("pO2_mmHg", "pO2", name)) %>%
  
  ggplot(aes(COPD_subclass, value, fill = COPD_subclass))+
  geom_dotplot(dotsize = 1.3, binaxis="y", stackdir="center") +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5)+
  
  labs(x = "", y = "Capillary partial pressure of oxygen \n(mmHg)")+
  stat_compare_means(aes(label = paste0(after_stat(p.signif))), size = 8,
                     symnum.args = symnum.args )+ 
  
  axis+
  theme(legend.position="none")+
  
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  scale_fill_manual(values=c("A" = "#661510", "B" = "#e34e45", "Donor"  = "gray75"))+
  scale_x_discrete(labels=c("A" = "Subgroup I", "B" = "Subgroup II"))+
  facet_wrap(~name, scales = "free_y", ncol = 2)

f5b2

patch <- f5b1+f5b2+ plot_layout(axes = "collect")
patch

save_plot(paste0("plots/", Sys.Date(), "_COPD_subclass_FACS_Serum-gas_exchange", ".png"), plot = patch, 
          base_height = 5.5, base_width = 8, dpi = 600)










