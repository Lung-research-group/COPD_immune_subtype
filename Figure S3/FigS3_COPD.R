# Figures S3 ---------------

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
facs.cc <- read.xlsx("iScience_input_data/hFACS_Donor_COPD_cc_full.xlsx") %>% 
  mutate(Unique_Sample_ID = as.factor(Unique_Sample_ID)) %>%
  mutate(Diagnosis = factor(Diagnosis, levels = c("Donor", "COPD")))

facs.cc.trans <- read.xlsx("iScience_input_data/hFACS_Donor_COPD_cc_LOG_x1c_full.xlsx") %>% 
  mutate(Unique_Sample_ID = as.factor(Unique_Sample_ID)) %>%
  mutate(Diagnosis = factor(Diagnosis, levels = c("Donor", "COPD")))

Parameter_labels <- read_excel("iScience_input_data/Parameter_labels.xlsx", 1)
parameter_labeller <- Parameter_labels$expression
names(parameter_labeller) <- Parameter_labels$Parameter_Name

# Graphs ------------------
## Redundant parental populations
double.d <- c("Lymphocytes", # parental pop
              "Macrophages", # parental pop
              "Monocytes", # parental pop
              "PMNL", # parental pop
              "Mast", # parental pop
              "Macs_CD14med", "Macs_CD14hi", # parental pop
              "DC_CD209neg_CD11cpos", 
              "DC_CD209pos_CD11cneg", 
              "DC_CD209pos_CD11cpos") 



# Fig S3A ------------------
# see https://r-graph-gallery.com/299-circular-stacked-barplot.html

bar.chart.pops <- c("Basophils", "Lymphocytes", "Macrophages", "Monocytes", "Neutrophils", "CD45pos_CD203neg_CD117pos", "CD45pos_CD203neg_CD117pos",
                    "DC_CD209neg_CD11cpos", "DC_CD209pos_CD11cneg", "DC_CD209pos_CD11cpos", 
                    "GRAN_CD193pos_CD16neg", "GRAN_CD193pos_CD16pos")

data <- facs.cc %>% pivot_longer(where(is.numeric)) %>% rename(individual = Unique_Sample_ID, group = Diagnosis) %>% filter(name %in% bar.chart.pops)  

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
levels(data$observation)



# Set a number of 'empty bar' to add at the end of each group
empty_bar=2
nObsType=nlevels(as.factor(data$observation))
to_add = data.frame( matrix(NA, empty_bar*nlevels(data$group)*nObsType, ncol(data)) )
colnames(to_add) = colnames(data)
to_add$group=rep(levels(data$group), each=empty_bar*nObsType )
data=rbind(data, to_add)
data=data %>% arrange(group, individual)
data$id=rep( seq(1, nrow(data)/nObsType) , each=nObsType)

# Get the name and the y position of each label
label_data= data %>% group_by(id, individual) %>% summarize(tot=sum(value))
number_of_bar=nrow(label_data)
angle= 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust<-ifelse( angle < -90, 1, 0)
label_data$angle<-ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data=data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data = base_data
grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start = grid_data$start - 1
grid_data=grid_data[-1,]


# Make the plot
p <- ggplot(data) +      
  
  # Add the stacked bar
  geom_bar(aes(x=as.factor(id), y=value, fill=observation), stat="identity", alpha=0.8) +
  #scale_fill_viridis(discrete=TRUE) +
  star_man.fill+
  labs(fill = "Cell type")+
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey60", alpha=1, size=0.5 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 2000, xend = start, yend = 2000), colour = "grey60", alpha=1, size=0.5 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 5000, xend = start, yend = 5000), colour = "grey60", alpha=1, size=0.5 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 10000, xend = start, yend = 10000), colour = "grey60", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 20000, xend = start, yend = 20000), colour = "grey60", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data$id),5), y = c(0, 2000, 5000, 10000, 20000), label = c("0", "2e4", "5e4", "1e5", "2e5") , 
           color="grey60", size=4 , angle=0, fontface="bold", hjust=1) +
  
  ylim(-12000,max(label_data$tot, na.rm=T)) +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() +
  
  # Add labels on top of each bar
  #geom_text(data=label_data, aes(x=id, y=tot+10, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=3, angle= label_data$angle, inherit.aes = FALSE ) +
  ## it appears that W224 has more than 100% of cells
  
  # Add base line information
  geom_segment(data=base_data[1,], aes(x = start, y = -5, xend = end, yend = -5), colour = "#791812", alpha=0.8, size=1.5 , inherit.aes = FALSE )  +
  geom_segment(data=base_data[2,], aes(x = start, y = -5, xend = end, yend = -5), colour = "#757b87", alpha=1, size=1.5 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -18, label=group), hjust=c(1,0), colour = "black", alpha=0.8, size=5, fontface="bold", inherit.aes = FALSE)+
  theme(legend.position = "right")
p
legend <- get_legend(p)
plot(legend)

save_plot(paste0("plots/", Sys.Date(), "_FACS_data_stacked_bar_cc", "_legend", ".png"), 
          plot = legend, base_height = 2, base_width = 1.5, dpi = 300)

p <- plot(p) +theme(legend.position = "none")


save_plot(filename = paste0("Plots/", Sys.Date(), '_FACS_data_stacked_bar_cc', '.png'), 
          plot = p, base_height = 7, base_width = 7, limitsize = FALSE)



###  quantification -------------------
## sig.df 
sig.df <- data %>%
  na.omit() %>%
  group_by(observation) %>% 
  summarise(pval = wilcox.test(value ~ group)$p.value) %>% 
  mutate(padj = p.adjust(pval, method = "BH", n = length(pval))) %>%
  mutate(pval = round(pval, 5), padj = round(padj, 5)) %>%
  mutate(log.padj = round(-log10(padj), 5)) %>%
  mutate(signif = case_when(padj <= 0.05 & padj > 0.01 ~"*", padj <= 0.01 & padj > 0.001 ~"**", padj <= 0.001 ~"***", .default = "ns"))

# Fig S3B ------------------
data %>% 
  filter(!is.na(value)) %>%
  rename(Diagnosis = group) %>%
  mutate(Diagnosis = factor(Diagnosis, levels = c("Donor", "COPD"))) %>%
  mutate(value = log10(value+1)) %>%
  
  ggplot(., aes(x = Diagnosis, y = value, fill = Diagnosis))+
  geom_violin(color="white", alpha=0.2, lwd= 1)+
  geom_dotplot(dotsize = 1, binaxis="y", stackdir="center") +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5)+
  labs(x = "",       
       y= expression("Cell number per "~mg^"-1"~tissue~(LOG~transformed)))+
  stat_compare_means(aes(label = paste0(after_stat(p.signif))), size = 8,
                     symnum.args = symnum.args )+ 
  man.fill+
  axis+
  theme(legend.position="none")+
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.2)))+
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
  theme(axis.text.x = element_text(size = 18))+
  facet_wrap(~observation, scales = "fixed", ncol = 5)

save_plot(paste0("plots/", Sys.Date(), "_FACS_data_pop_overview_cc", ".png"), 
          plot = ggplot2::last_plot(), base_height = 5, base_width = 12, dpi = 600)


cells = c("CD3", "CD4", "CD19", "Mono_int","Neutrophils", "CD8")
cell.order.w = cells

###  sig.df
sig.df <- facs.cc.trans %>% 
  #select(Unique_Sample_ID, Diagnosis, all_of(cells)) %>% 
  pivot_longer(where(is.numeric)) %>%
  na.omit(value) %>%
  group_by(name) %>%
  summarise(pval = wilcox.test(value ~ Diagnosis)$p.value) %>% 
  mutate(padj = p.adjust(pval, method = "BH", n = length(pval))) %>%
  mutate(pval = round(pval, 5), padj = round(padj, 5)) %>%
  mutate(log.padj = round(-log10(padj), 5)) %>%
  mutate(signif = case_when(padj <= 0.05 & padj > 0.01 ~"*", padj <= 0.01 & padj > 0.001 ~"**", padj <= 0.001 ~"***", .default = "ns"))

# Fig S3H ------------------
facs.cc.trans %>% 
  select(Unique_Sample_ID, Diagnosis, all_of(cells)) %>% 
  pivot_longer(all_of(cells)) %>% 
  mutate(name = factor(name, levels = cell.order.w)) %>% 
  
  ggplot(., aes(x = Diagnosis, y = value, fill = Diagnosis))+
  geom_violin(color="white", alpha=0.2, lwd= 1)+
  geom_dotplot(dotsize = 2, binaxis="y", stackdir="center") +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5)+
  labs(x = "",       
       y= expression("Cell number per "~mg^"-1"~tissue~(Log[10]~transformed)))+
  #stat_compare_means(aes(label = paste0(after_stat(p.signif))), size = 8, symnum.args = symnum.args ) + 
  man.fill+
  man.col+
  axis+
  theme(legend.position="none")+
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.25)))+
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
  facet_wrap(~name, scales = "fixed", labeller=as_labeller(parameter_labeller, label_parsed), ncol = 6)

save_plot(paste0("plots/", Sys.Date(), "_FACS-v10_RFpops_cc_w_reduced", ".png"), 
          plot = ggplot2::last_plot(), base_height = 8, base_width = 10, dpi = 600)


## all remaining pops --------------
facs.cc.trans %>% 
  select(Unique_Sample_ID, Diagnosis, !all_of(double.d)) %>% #removes redundant parental populations
  select(!all_of(cells)) %>% #removes cells plotted in Fig 1H
  pivot_longer(where(is.numeric)) %>% 
  mutate(name = gsub("NKcells", "NK.cells", name)) %>%
  mutate(Parameter_Name = name) %>%
  left_join(select(Parameter_labels, Parameter_Name, expression)) %>% 
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
