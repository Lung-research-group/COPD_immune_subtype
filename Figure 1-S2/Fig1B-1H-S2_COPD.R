# Figures 1 and S2 ---------------

# Setup ----------------
## Load the required libraries --------------
library(tidyverse)
library(ggpubr)
library(openxlsx)
library(ggpubr)
library(cowplot)
require(ggrepel)
library(patchwork)
library(readxl)
library(tidyr)

dir.create("Plots")

## Formatting ---------------
theme_set(theme_pubr(base_size=16, border = T))

man.fill <-   scale_fill_manual(values = c("Donor" = "#757b87", "COPD" = "#791812"))
man.col <-   scale_color_manual(values = c("Donor" = "#757b87", "COPD" = "#791812"))

star_man.fill <- scale_fill_manual(
  values=c("DC" = "#80C490", 
           "Macrophages" = "#64b2ce", 
           "Monocytes" = "#80dcc5", 
           "Lymphocytes"  = "#fc5361",
           "PMNL" = "#5e738f"
  ))

axis <- theme(axis.line.x = element_line(color="black", size = 0.75),
              axis.line.y = element_line(color="black", size = 0.75), 
              plot.title = element_text(size = rel(1), colour="black", element_text(hjust = 0.5)),
              
              axis.title.x = element_text(margin=margin(0,10,0,0)), 
              axis.title.y = element_text(margin=margin(0,10,0,0)), 
              axis.text = element_text(colour="black"), 
              
              panel.grid = element_blank(),
              plot.caption = element_text(size = rel(0.6)))


symnum.args <- list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), 
                    symbols = c("***", "**", "*", "ns"))


# File import --------------
facs.per <- read.xlsx("iScience_input_data/hFACS_Donor_COPD_per_full.xlsx") %>% mutate(Unique_Sample_ID = as.factor(Unique_Sample_ID))
facs.per.trans <- read.xlsx("iScience_input_data/hFACS_Donor_COPD_per_LOG_x1c_full.xlsx") %>% mutate(Unique_Sample_ID = as.factor(Unique_Sample_ID))

Parameter_labels <- read_excel("iScience_input_data/Parameter_labels.xlsx", 1)
parameter_labeller <- Parameter_labels$expression
names(parameter_labeller) <- Parameter_labels$Parameter_Name

# Graphs ------------------
## Redundant parental populations ---------
double.d <- c("Lymphocytes", # parental pop
              "Macrophages", # parental pop
              "Monocytes", # parental pop
              "PMNL", # parental pop
              "Mast", # parental pop
              "Macs_CD14med", "Macs_CD14hi", # parental pop
              "DC_CD209neg_CD11cpos", 
              "DC_CD209pos_CD11cneg", 
              "DC_CD209pos_CD11cpos") 
              
# Figure 1B --------------
## Stacked areaplot ----------------
bar.chart.pops <- c("Basophils", "Lymphocytes", "Macrophages", "Monocytes", "Neutrophils", 
                    "CD45pos_CD203neg_CD117pos", "CD45pos_CD203neg_CD117pos",
                    "DC_CD209neg_CD11cpos", "DC_CD209pos_CD11cneg", "DC_CD209pos_CD11cpos", 
                    "GRAN_CD193pos_CD16neg", "GRAN_CD193pos_CD16pos")

data <- facs.per %>% mutate(Lymphocytes = Lymphocytes-Basophils) %>% #Basophils removed from total lymph counts
  pivot_longer(where(is.numeric)) %>% rename(individual = Unique_Sample_ID, group = Diagnosis) %>% 
  filter(name %in% bar.chart.pops)  

data <- data %>% mutate(Parent = case_when(grepl("Basophils", name) ~ "PMNL",
                                           grepl("Lymphocytes", name) ~ "Lymphocytes",
                                           grepl("Macrophages", name) ~ "Macrophages",
                                           grepl("Monocytes", name) ~ "Monocytes",
                                           grepl("Neutrophils", name) ~ "PMNL",
                                           grepl("CD45pos_CD203", name) ~ "PMNL",
                                           grepl("DC_CD209", name) ~ "DC",
                                           grepl("GRAN_CD193", name) ~ "PMNL"))

data <- data %>% group_by(individual, group, Parent) %>% summarize(value=sum(value, na.rm = TRUE)) %>% rename(observation = Parent)

data <- data %>% mutate(observation = factor(observation, levels = c("DC", "Macrophages", "Monocytes", "Lymphocytes", "PMNL")))

data <- data %>% group_by(individual) %>% mutate(prop = prop.table(value)) %>% mutate(value.orig = value,
                                                                                      value = prop*100)

data %>% 
  group_by(group, observation) %>% 
  mutate(numbering =  seq_along(individual)) %>%
  mutate(group = factor(group, levels = c("Donor", "COPD"))) %>%
  mutate(., x = ifelse(is.na(prop), 0, prop)) %>%
  mutate(group.id = paste0(group, ".", individual)) %>%
  
  ggplot(., aes(x=numbering, y=prop, fill=observation)) + 
  
  geom_area(alpha=0.6 , size=.5, colour="white") +
  star_man.fill+
  labs(y="Proportion")+
  scale_x_continuous(breaks = 1:23)+
  scale_y_continuous(expand = expansion(0), labels = scales::percent_format()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  axis+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        panel.grid.minor.x = element_line(colour = "grey", size = 0.5),
        legend.position = "none")+

  facet_wrap(~group, scales = "free_x")


save_plot(paste0("plots/", Sys.Date(), "_FACS_data_pop_overview0", ".png"), 
          plot = ggplot2::last_plot(), base_height = 5, base_width = 8, dpi = 600)



## sig.df ---------
sig.df <- data %>%
  na.omit() %>%
  group_by(observation) %>% 
  summarise(pval = wilcox.test(value ~ group)$p.value) %>% 
  mutate(padj = p.adjust(pval, method = "BH", n = length(pval))) %>%
  mutate(pval = round(pval, 5), padj = round(padj, 5)) %>%
  mutate(log.padj = round(-log10(padj), 5)) %>%
  mutate(signif = case_when(padj <= 0.05 & padj > 0.01 ~"*", padj <= 0.01 & padj > 0.001 ~"**", padj <= 0.001 ~"***", .default = "ns"))


# Figure S2A -------------
## Relative changes in global distribution -----------------------
data %>% 
  mutate(group = factor(group, levels = c("Donor", "COPD"))) %>%
  filter(!is.na(value)) %>%
  rename(Diagnosis = group) %>%
  mutate(value = log10(value+0.0001)) %>%
  
  ggplot(., aes(x = Diagnosis, y = value, fill = Diagnosis))+ 
  geom_dotplot(dotsize = 1.5, binaxis="y", stackdir="center") +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5)+
  labs(x = "", y= expression("%"*CD45^"+"~cells~(LOG~transformed)))+
  stat_compare_means(aes(label = paste0(after_stat(p.signif))), size = 8,
                     symnum.args = symnum.args )+ 
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.25)))+
  man.fill+
  axis+
  theme(legend.position="none",
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), 
        axis.text.x = element_text(size = 18))+
  facet_wrap(~observation, scales = "fixed", ncol = 5)

save_plot(paste0("plots/", Sys.Date(), "_FACS_data_pop_overview", ".png"), 
          plot = ggplot2::last_plot(), base_height = 5, base_width = 12, dpi = 600)




# Figure 1C-F --------------
## see the following scripts -------------
##R_Fig1C_PCA_biplot.r
##R_Fig1D_OPLSDA_Scores_plot.r
##R_Fig1E_F_G_RandomForest_Vardepth_MDS_plots.r

# Figure 1H  --------
## plot 7 most important cells ---------
cells = c("Mono_clas", "CD19", "CD3", "Mono_int", "CD4", "CD8", "Neutrophils")
cell.order.w = cells

facs.per.trans %>% 
  pivot_longer(all_of(cells)) %>% 
  mutate(name = factor(name, levels = cell.order.w)) %>% 
  
  ggplot(., aes(x = Diagnosis, y = value, fill = Diagnosis))+
  geom_violin(color="white", alpha=0.2, lwd= 1)+
  geom_dotplot(dotsize = 1.5, binaxis="y", stackdir="center") +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5)+
  labs(x = "", y= expression("%"*CD45^"+"~cells~(LOG~transformed)))+
  stat_compare_means(aes(label = paste0(after_stat(p.signif))), size = 8,
                     symnum.args = symnum.args )+ 
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.25)))+
  man.fill+
  axis+
  theme(legend.position="none",
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        axis.text.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)) + 
  facet_wrap(~name, scales = "fixed", labeller=labeller(name = c("Mono_clas" = "Monocytes \nclassical",  "CD19" = "B cells",       
                                                                 "CD3" = "T cells", "Mono_int" ="Monocytes \nintermediate" ,   
                                                                 "CD4" = "T helper", "CD8" = "T cytotoxic", "Neutrophils" ="Neutrophils")), ncol = 7)+
  theme(strip.text = element_text(size = 18))

save_plot(paste0("plots/", Sys.Date(), "_FACS_RFpops_top7", ".png"), 
          plot = ggplot2::last_plot(), base_height = 5, base_width = 14, dpi = 600)


# Figure S2B ------------
## NLR per ---------
facs.per %>% 
  mutate(NLR = PMNL/Lymphocytes) %>%
  select(Unique_Sample_ID, Diagnosis, NLR) %>% 
  pivot_longer(where(is.numeric)) %>% 

  ggplot(., aes(x = Diagnosis, y = value, fill = Diagnosis))+
  geom_violin(color="white", alpha=0.4, lwd= 1)+
  geom_dotplot(dotsize = 1.5, binaxis="y", stackdir="center") +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5)+
  labs(x = "",  y = "Neutrophil to lymphocyte ratio")+ 
  stat_compare_means(aes(label = paste0(after_stat(p.signif))), size = 8, symnum.args = symnum.args ) + 
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.25)))+
  man.fill+
  axis+
  theme(legend.position="none", 
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
  facet_wrap(~name, scales = "fixed")

save_plot(paste0("plots/", Sys.Date(), "_FACS_NLR_CD45", ".png"), 
          plot = ggplot2::last_plot(), base_height = 5, base_width = 3.25, dpi = 600)


# Figure S2C ----------------
## All remaining FACS pops -----
facs.per.trans %>% 
  select(Unique_Sample_ID, Diagnosis, !all_of(double.d)) %>% #removes redundant parental populations
  select(!all_of(cells)) %>% #removes cells plotted in Fig 1H
  pivot_longer(where(is.numeric)) %>% 
  mutate(name = gsub("NKcells", "NK.cells", name)) %>%
  mutate(Parameter_Name = name) %>%
  left_join(select(Parameter_labels, Parameter_Name, expression)) %>% #Parameter_labels file needed
  mutate(Diagnosis = factor(Diagnosis, levels = c("Control", "COPD"))) %>%
  
  ggplot(aes(value, expression, fill = Diagnosis, colour = Diagnosis)) +
  geom_violin(alpha=0.4, lwd= 0.5, position = position_dodge(width = 0.75), scale = "count", adjust = 2)+
  geom_point(position = position_dodge(width = 0.75), colour = "black", size = 3)+
  geom_point(position = position_dodge(width = 0.75), size = 2)+
  stat_summary(fun = median, fun.min = median, fun.max = median, position = position_dodge(width = 0.75),
               geom = "crossbar", width = 0.5, colour = "black")+
  man.fill+
  man.col+
  labs(y = "",   subtitle = expression("Percentage CD45"^"+"~"cells"), x= expression("%"*CD45^"+"~cells~(LOG~transformed)))+
  scale_y_discrete(labels = scales::label_parse(), limits=rev)+
  theme(legend.position = "none")

save_plot(paste0("plots/", Sys.Date(), "_FACS_AllPops_per_remaining", ".png"), 
          plot = ggplot2::last_plot(), base_height = 9, base_width = 8, dpi = 600)


## sig.df ---------
sig.df <- facs.per.trans %>% 
  select(Unique_Sample_ID, Diagnosis, !all_of(double.d)) %>% #removes redundant parental populations
  select(!all_of(cells)) %>% #removes cells plotted in Fig 1H
  pivot_longer(where(is.numeric)) %>% 
  mutate(name = gsub("NKcells", "NK.cells", name)) %>%
  na.omit() %>%
  group_by(name) %>% 
  summarise(pval = wilcox.test(value ~ Diagnosis)$p.value) %>% 
  mutate(padj = p.adjust(pval, method = "BH", n = length(pval))) %>%
  mutate(pval = round(pval, 5), padj = round(padj, 5)) %>%
  mutate(log.padj = round(-log10(padj), 5)) %>%
  mutate(signif = case_when(padj <= 0.05 & padj > 0.01 ~"*", padj <= 0.01 & padj > 0.001 ~"**", padj <= 0.001 ~"***", .default = "ns"))




# Figure S2D ----------------
## CD4/CD8 per --------------------
facs.per %>% 
  mutate(CD4.CD8 = CD4/CD8, NLR = PMNL/Lymphocytes) %>%
  select(Unique_Sample_ID, Diagnosis, CD4.CD8) %>% 
  pivot_longer(where(is.numeric)) %>% 
  mutate(name = gsub("CD4.CD8", "%CD45 cells", name)) %>%

  ggplot(., aes(x = Diagnosis, y = value, fill = Diagnosis))+
  geom_violin(color="white", alpha=0.4, lwd= 1)+
  geom_dotplot(dotsize = 1.5, binaxis="y", stackdir="center") +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5)+
  labs(x = "",  subtitle =  "", y = "CD4/CD8 ratio")+ 
  stat_compare_means(aes(label = paste0(after_stat(p.signif))), size = 8, symnum.args = symnum.args ) + 
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.25)))+
  man.fill+
  man.col+
  axis+
  theme(legend.position="none", 
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
  facet_wrap(~name, scales = "fixed")

save_plot(paste0("plots/", Sys.Date(), "_FACS_CD4-CD8ratio", ".png"), 
          plot = ggplot2::last_plot(), base_height = 6, base_width = 3.5, dpi = 600)



