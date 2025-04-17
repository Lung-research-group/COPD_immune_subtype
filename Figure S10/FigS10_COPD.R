# Figures S10 ---------------

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


# import data ---------
patho <- read.xlsx("iScience_input_data/patho.xlsx") %>% mutate(Unique_Sample_ID = as.factor(Unique_Sample_ID))

# make column factors
cols <- names(patho[3:ncol(patho)])
patho[cols] <- lapply(patho[cols], factor) 



##Centracinar.Empysema ------------------------
m <- MASS::polr(Centracinar.Empysema ~ COPD_subclass, data = patho, Hess = TRUE, method = "logistic")

## view a summary of the model
summary(m)
## store table
(ctable <- coef(summary(m)))

## calculate and store p values
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2

## combined table
(ctable <- cbind(ctable, `p value` = round(p, 3)))


## Respiratory bronchiolitis --------
fisher.test(table(patho$COPD_subclass, patho$Respiratory.bronchiolitis)) # p-value = 1

## Intraluminal inflammatory exudate
fisher.test(table(patho$COPD_subclass, patho$Intraluminal.inflammatory.exudate)) # p-value = 0.4947


## Intramural chronic inflammation
m <- MASS::polr(Intramural.chronic.inflammation ~ COPD_subclass, data = patho, Hess = TRUE, method = "logistic")
summary(m)
(ctable <- coef(summary(m)))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
(ctable <- cbind(ctable, `p value` = round(p, 3)))


## Lymphoid.follicles -----------
m <- MASS::polr(Lymphoid.follicles ~ COPD_subclass, data = patho, Hess = TRUE, method = "logistic")
summary(m)
(ctable <- coef(summary(m)))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
(ctable <- cbind(ctable, `p value` = round(p, 3)))


## BALT  ------------------------
m <- MASS::polr(BALT.hyperlasia ~ COPD_subclass, data = patho, Hess = TRUE, method = "logistic")
summary(m)
(ctable <- coef(summary(m)))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
(ctable <- cbind(ctable, `p value` = round(p, 3)))


## Smooth muscle hypetrophy	 ---------------
fisher.test(table(patho$COPD_subclass, patho$Smooth.muscle.hypetroph)) # p-value = 0.1287

## Goblet.cell.metaplasia -------------
m <- MASS::polr(Goblet.cell.metaplasia ~ COPD_subclass, data = patho, Hess = TRUE, method = "logistic")
summary(m)
(ctable <- coef(summary(m)))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
(ctable <- cbind(ctable, `p value` = round(p, 3)))



## graph --------------------------
patho %>%
  select(COPD_subclass, COPD_subclass, Centracinar.Empysema:Goblet.cell.metaplasia, 
         -Squamous.metaplasia, -BALT.hyperlasia, -Subepithelial.fibrosis, -Intraluminal.inflammatory.exudate, -Adventitial.fibrosis) %>%
  pivot_longer(Centracinar.Empysema:Goblet.cell.metaplasia) %>%
  mutate(status = case_when(value == 0 ~ "Absent",
                            value == 1 ~ "Minor",
                            value == 2 ~"Medium",
                            value == 3 ~ "Severe")) %>%
  mutate(status = factor(status, levels = c("Severe", "Medium", "Minor", "Absent"))) %>%
  na.omit %>%
  
  mutate(name = gsub("Centracinar.Empysema", "Centracinar emphysema", name),
         name = gsub("Lymphoid.follicles", "Lymphoid follicles", name),
         name = gsub("Intramural.chronic.inflammation", "Intramural inflammation", name),
         name = gsub("Respiratory.bronchiolitis","Respiratory bronchiolitis" , name),
         name = gsub("Goblet.cell.metaplasia", "Goblet cell metaplasia", name),
         name = gsub("Smooth.muscle.hypetrophy", "Smooth muscle hypertrophy", name)) %>%
  
  mutate(name = factor(name, levels = c("Centracinar emphysema", "Lymphoid follicles", 
                                        "Intramural inflammation", "Respiratory bronchiolitis", 
                                        "Goblet cell metaplasia" , "Smooth muscle hypertrophy" ))) %>%
  
  group_by(COPD_subclass, name) %>%
  count(name, status) %>%
  mutate(prop = prop.table(n)) %>%
  
  ggplot(., 
         aes(x=COPD_subclass, 
             y=prop, 
             fill=status)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("Severe" = "#b30000", "Medium" = "#e34a33", "Minor"  = "#fc8d59","Absent" = "#fdcc8a"))+
  labs(x = "", 
       title = " ")+
  scale_y_continuous(expand = expansion(mult = c(NA, 0.1)), labels = scales::percent, breaks = c(0, 0.25, .50, 0.75, 1.00))+
  scale_x_discrete(labels=c("A" = "Subgroup \nI", "B" = "Subgroup \nII"))+
  labs(y = "Relative proportion", x = "") +
  theme(legend.position = "right")+
  facet_wrap(~name, scales = "fixed", ncol = 6)

save_plot(paste0("plots/", Sys.Date(), "_Pathology_ConSubclass", ".png"), 
          plot = ggplot2::last_plot(), base_height = 5.5, base_width = 16.5, dpi = 600)






