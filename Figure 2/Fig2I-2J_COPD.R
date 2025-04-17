#libraries ====
library(readxl)
library(ggpubr)
library(ggplot2)
library(igraph)
library(effsize)
library(reshape2)
library(dplyr)
library(stringr)
library(ggforce)
library(ggrepel)
library(ggpp)
library(ggforce)

#rev 01 approach (2) only COPD data --------------------------
dir <- paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/")
final_data <- read.csv(file =paste0(dir, "input data/final_data_revSC_COPDvsDonor.csv"))
final_data$Consensus_subclass <- NULL # remove the subclass diagnosis
elisa_new <- as.data.frame(final_data) # convert to data frame
elisa_new <- elisa_new[which(elisa_new$Diagnosis %in% c("COPD")),] #ONLY TAKING PATIENTS WITH IMMUNE PROFILE
clinical <- c("Age","BMI","Smoking_py","CRP",
              "FEV1_percent","FVC_percent","FEV1_FVC_percent",          
              "RV_percent","mPAP","DLCOcSB_percent","pO2_mmHg",                  
              "pCO2_mmHg")
sig_parameters <- c("CD19","CD3","CD4","CD8",
                    "Neutrophils" , "Macs_CD14hi_CD1aposHLAneg",
                    "Mono_clas","Mono_int", "CCL5_L","CXCL10_L",
                    "CXCL9_L","IL8_L",
                    "CXCL1_S","CCL17_S","CCL4_S","IL1b_S",
                    "CCL2_S","IL10_S","CXCL5_S","CCL5_S",
                    "IL6_S","CXCL9_S","CXCL10_S",
                    "TNFb_S","IFN_la2_3_S" ,clinical)

up_parameter <- c("CD19","CD3","CD4","CD8","CXCL10_L","CXCL1_S","CCL17_S",
                  "CXCL9_L", "CCL5_L","CCL4_S","IL1b_S", "CCL2_S","IL10_S",
                  "CXCL5_S","CCL5_S","IL6_S")
down_parameter <- c("Neutrophils",
                    "Macs_CD14hi_CD1aposHLAneg",
                    "Mono_clas","Mono_int", "IL8_L",
                    "CXCL9_S","CXCL10_S","TNFb_S","IFN_la2_3_S")
elisa_new <- elisa_new[which(colnames(elisa_new) %in% sig_parameters)]


temp <- colnames(elisa_new)# we want to convert to numeric all columns in the data
for (i in temp[!temp == c("Diagnosis")]){
  elisa_new[,i] <- as.numeric(elisa_new[,i])
}

#prepare the parameter excel file. 
parameter <- read_xlsx(paste0(dir, "input data/parameter_label.xlsx")) # it should be easily found inside the data folder. Here inside the parameter label and the corresponding labels

#calculate fold change (COPD vs DOnor) for each parameter
clinical <- c("Age","BMI","Smoking_py","CRP",
              "FEV1_percent","FVC_percent","FEV1_FVC_percent",          
              "RV_percent","mPAP","DLCOcSB_percent","pO2_mmHg",                  
              "pCO2_mmHg")
cells <- c( "Basophils", "CD19","CD3",                       
            "CD4","CD8","DC_CD209neg_CD11cposCd1a",
            "DC_CD209pos_CD11cposCD1a","gdTCR","Macs_CD14hi_CD1anegHLApos",
            "Mono_clas", "Mono_int","Mono_non",
            "Macs_CD14hi_CD1aposHLAneg",  "Macs_CD14hi_CD1aposHLApos", 
            "Macs_CD14med_CD1anegHLApos","Macs_CD14med_CD1aposHLAneg", 
            "Macs_CD14med_CD1aposHLApos","Monocytes","Neutrophils","Mast_CD203neg","Mast_CD203pos",
            "NKcells", "NKT","pDC","Tregs") 
cytokine <- colnames(final_data)[-which(colnames(final_data) %in% c(clinical, cells, "Unique_Sample_ID","Diagnosis", "Consensus_subclass"))]

colnames(elisa_new)
final_colnames_elisa <- colnames(elisa_new)
estrogen <- as.matrix(elisa_new[, final_colnames_elisa]) # we filter the data only  keeping the significant parameters

#get the corr matrix and p val matrix
library("Hmisc")
corr_df <- rcorr(estrogen, type="pearson")
dt <- corr_df$r
dp <- corr_df$P

g <- graph.adjacency( # Create a graph adjacency based on correlation distances between parameters in  pairwise fashion.
  as.matrix(as.dist(dt)), # so basically those which correlated highly will be positioned closer in the graph
  mode="undirected",
  weighted=TRUE, # we want to weigh the edge (weight here is the correlation coefficient)
  diag=FALSE
)


g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE)# Simplfy the adjacency object


g <- delete_edges(g, E(g)[which(abs(E(g)$weight)<0.5)])# Remove edges below absolute Pearson correlation 0.5 and with NA value
g <- delete_edges(g, E(g)[!complete.cases(E(g)$weight)])


V(g)$name <- V(g)$name # Assign names to the graph vertices (optional). vertices here is the parameters name

y <- get.data.frame(g) # get the igraph object with weight of each connection
y$change <- case_when(y$weight > 0 ~ "positive corr", # annotate the weight or the correlation coeficient
                      y$weight < 0 ~ "negative corr")
y$corr <- paste0(y$from, "_", y$to) # label the connecting parameters
E(g)$weight <- abs(E(g)$weight) # add the weight into edges information. Edges here is the connecting parameters

#create metadata
set.seed(1234)
fr.all <- layout.fruchterman.reingold(g) # we apply fruchterman reingold to draw the graph
fr.all.df <- as.data.frame(fr.all) # we get the coordinate x and y as V1 and V2 columns
fr.all.df$species <- final_colnames_elisa #we attach parameters name into coordinate values

#create correlation data frame to be saved as supplementary
x <- get.data.frame(g)
x$corr <- paste0(x$from, "_", x$to)
x$change <- y$change[match(x$corr, y$corr)] #attach the information from y table above, the information about positiive and negative corr
x$from.x <- fr.all.df$V1[match(x$from, fr.all.df$species)]*10  #we took V1 parameters from metadata above, and put them here in table as from.x column;match with paramater in  column "from" to  column "species" in metadata for coordinate x (V1)
x$from.y <- fr.all.df$V2[match(x$from, fr.all.df$species)]*10 # match with paramater in  column "from" to  column "species" in metadata for coordinate y (V2)

x$to.x <- fr.all.df$V1[match(x$to, fr.all.df$species)]*10  #match with paramater in  column "to" to  column "species" in metadata for coordinate y (V1)
x$to.y <- fr.all.df$V2[match(x$to, fr.all.df$species)]*10 #match with paramater in  column "to" to  column "species" in metadata for coordinate y (V2)

fr.all.df$V1 <- fr.all.df$V1*10 # copy the coordinate and increase the coordinate as we did previously to other connections
fr.all.df$V2 <- fr.all.df$V2*10

x$change <- relevel(as.factor(x$change), ref="positive corr") # we relevel the factor of correlation label


temp <- as.data.frame(table(x$from)) # we calculate connections "from". we have two connections: from and to
colnames(temp)[2] <- "num_of_connection_Var1" # we label the name of connections
dat <- y[which(y$from %in% temp$Var1),] # filter data, only including the confirmed connections from both x and y (just as double check)
dat <- merge(dat, temp, by.x="from", by.y="Var1") # attach the rest information


#cluster the parameters using the fast greedy algorithm
fc <- cluster_fast_greedy(g)# we take the igraph object
cluster <- cut_at(fc, 5)
mem <-membership(fc)# assign clusters
z <- names(mem)#take all parameters name
y <- as.numeric(cluster)# take all cluster number
h <-data.frame(name =z, community =y)# make data frame of cluster name and parameters name
fr.all.df <- merge(fr.all.df, h, by.x="species", by.y="name")# attach the cluster name and parameters name to metadata we created above

a <- as.data.frame(table(dat$to)) # we take parameters name in column "to" from correlation data frame, and calculate how many times they appear indicating how many connections they have as from/origin
b <- as.data.frame(table(dat$from)) # we take parameters name in column "from" from correlation data frame, and calculate how many times they appear indicating how many connections they have as to/destination
hh <- full_join(a,b, by="Var1")  # we join the information of above a and b
hh[is.na(hh)] <- 0 # we assign any NAs to 0
hh$total <- hh$Freq.x+hh$Freq.y # we sum the total number of connections
fr.all.df <- merge(fr.all.df, hh, by.x="species", by.y="Var1") # and attach this information to metadata

colnames(parameter)[1] <- "species"
fr.all.df <- merge(fr.all.df, parameter, by="species")
fr.all.df$label <- ordered(fr.all.df$label, levels=c("clinical and lung function","circulation","lung")) # we also reorder the level of origin label of parameters


dp_sig <- melt(dp)
colnames(dp_sig) <- c("to", "from", "pval")
dp_sig$corr<- paste0(dp_sig$from, "_", dp_sig$to)
corr_var <- dp_sig$corr[which(dp_sig$pval <= 0.05)]

dat <- dat[which(dat$corr %in% corr_var),]
pmrt <- unique(c(dat$from,dat$to))


#Plotting
#with cluster annotations
theme_journal <- theme(
  axis.text.x = element_blank(),  # remove x-axis text
  axis.text.y = element_blank(), # remove y-axis text
  axis.ticks = element_blank(),  # remove axis ticks
  axis.title.x = element_blank(), # remove x-axis labels
  axis.title.y = element_blank(), # remove y-axis labels
  panel.background =  element_rect(fill='transparent'), 
  plot.background = element_rect(fill='transparent', color=NA),
  panel.border =element_blank(), 
  legend.position = "none",
  plot.title = element_blank(),
  text = element_text(size=16),
  panel.grid.major = element_blank(),  #remove major-grid labels
  panel.grid.minor = element_blank())  #remove minor-grid labels

x <- x[which(x$corr %in% corr_var),]
fr.all.df <- fr.all.df[which(fr.all.df$species %in% pmrt),]
fr.all.df$direction <- NA
all <- c(up_parameter, down_parameter)
fr.all.df$direction[which(fr.all.df$species %in% up_parameter)] <- "up"
fr.all.df$direction[which(fr.all.df$species %in% down_parameter)] <- "down"
fr.all.df$direction[-which(fr.all.df$species %in% all)] <- "none"
plot1 <- ggplot() +
  geom_link(data=x,aes(x=from.x,xend = to.x, y=from.y,yend = to.y, size=weight, linetype=change, alpha=weight),color="black") + # add line type
  geom_point(data=fr.all.df,aes(x=V1,y=V2,fill=as.factor(label),
                                color=as.factor(label), shape=direction),show_guide=T,
             size=4)+ # add colour scaling for group membership
  scale_shape_manual(values=c(24,25,21), breaks = c("up","down","none"))+
  theme_set(theme_pubr(base_size=20, border = T))+
  scale_alpha(range = c(0.4, 0.8))+
  scale_size(range = c(0,1))+
  scale_x_continuous(expand=c(0,4))+  # expand the x limits 
  scale_y_continuous(expand=c(0,4))+ # expand the y limits
  theme_bw()+  # use the ggplot black and white theme
  theme_journal+
  scale_fill_manual(values=c("red", "royalblue4", "goldenrod"), 
                    breaks=c("circulation","lung","clinical and lung function"))+
  scale_color_manual(values=c("red", "royalblue4", "goldenrod"), 
                     breaks=c("circulation","lung","clinical and lung function"))+
  geom_text_repel(data=fr.all.df,aes(x=V1,y=V2,label=expression), max.overlaps = Inf,parse=T );plot1

plot1
ggsave(filename = paste0( "rev1_copd_sigparameters.png"), 
       plot=plot1, width = 8.5, height = 8.5, units = c("in"), dpi = 300, scale = 0.8)

plot1 <- ggplot() +
  geom_link(data=x,aes(x=from.x,xend = to.x, y=from.y,yend = to.y, 
                       alpha=weight,size=weight, linetype=change),color="black") + # add line type
  geom_point(data=fr.all.df,aes(x=V1,y=V2,fill=as.factor(community), shape=direction),show_guide=T, type=21, size=4)+ # add colour scaling for group membership
  scale_shape_manual(values=c(24,25,21), breaks = c("up","down","none"))+
  theme_set(theme_pubr(base_size=20, border = T))+
  scale_alpha(range = c(0.4, 0.8))+
  scale_size(range = c(0,1))+
  scale_x_continuous(expand=c(0,4))+  # expand the x limits 
  scale_y_continuous(expand=c(0,4))+ # expand the y limits
  theme_bw()+  # use the ggplot black and white theme
  theme_journal+
  scale_fill_manual(values=c("darkgreen", "salmon"), 
                    breaks=c("1","2"))+
  geom_text_repel(data=fr.all.df,aes(x=V1,y=V2,label=expression),parse=T);plot1



#visualise the corr plot--------
corr_var_to <- as.character(dp_sig$to[which(dp_sig$pval <= 0.05)])
corr_var_from <- as.character(dp_sig$from[which(dp_sig$pval <= 0.05)])
main_effects <- unique(c(corr_var_to, corr_var_from))
estrogen <- as.matrix(elisa_new[ ,c(main_effects)])

#1)
parameters <- c("DLCOcSB_percent","CXCL9_S", 
                "CXCL10_L","FEV1_FVC_percent",
                "CXCL9_L","CD8", 
                "pCO2_mmHg", "CD8",
                "FEV1_percent", "Mono_int",
                "TNFb_S","CD8")

mtcars <- as.data.frame(estrogen[,which(colnames(estrogen) %in% parameters)])
df <- mtcars[, c("DLCOcSB_percent","CXCL9_S")]
df <- na.omit(df)
df$DLCOcSB_percent <- as.numeric(df$DLCOcSB_percent)
df$CXCL9_S <- as.numeric(df$CXCL9_S)

a <- ggplot(df, aes(x=CXCL9_S, y=DLCOcSB_percent)) + 
  geom_point(color="#791812", size=4)+
  geom_smooth(method="lm")+
  ylab("DLCOcSB%")+xlab("sCXCL9")+
  theme(strip.text.x = element_text(size = 18),
        axis.text.x = element_text(size = 18, angle = 0, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        text = element_text(size = 18))+ggtitle("")+
  theme_bw(base_size = 18)+
  theme(text=element_text(size=18,  family="Arial"))#+expand_limits(y = 0)
a_plot <- a+theme(legend.position = "none")
a_plot
ggsave(filename = paste0("correlation_COPD_DLCOcSB_sCXCL9.png"), 
       plot=a_plot, width = 4, height = 4, units = c("in"), dpi = 300, scale = 1)

#2)
parameters <- c("DLCOcSB_percent","CXCL9_S", 
                "CXCL10_L","FEV1_FVC_percent",
                "CXCL9_L","CD8", 
                "pCO2_mmHg", "CD8",
                "FEV1_percent", "Mono_int",
                "TNFb_S","CD8")

mtcars <- as.data.frame(estrogen[,which(colnames(estrogen) %in% parameters)])
df <- mtcars[, c("FEV1_FVC_percent","CXCL10_L")]
df <- na.omit(df)
df$FEV1_FVC_percent <- as.numeric(df$FEV1_FVC_percent)
df$CXCL10_L <- as.numeric(df$CXCL10_L)

a <- ggplot(df, aes(x=CXCL10_L, y=FEV1_FVC_percent)) + 
  geom_point(color="#791812", size=4)+
  geom_smooth(method="lm")+ylab("FEV1/FVC%")+xlab("CXCL10")+
  theme(strip.text.x = element_text(size = 18),
        axis.text.x = element_text(size = 18, angle = 0, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        text = element_text(size = 18))+ggtitle("")+
  theme_bw(base_size = 18)+
  theme(text=element_text(size=18,  family="Arial"))#+expand_limits(y = 0)
b_plot <- a+theme(legend.position = "none")
b_plot
ggsave(filename = paste0("correlation_COPD_FEVFVC_CXCL10.png"), 
       plot=a_plot, width = 4, height = 4, units = c("in"), dpi = 300, scale = 1)


#3)
parameters <- c("DLCOcSB_percent","CXCL9_S", 
                "CXCL10_L","FEV1_FVC_percent",
                "CXCL9_L","CD8", 
                "pCO2_mmHg", "CD8",
                "FEV1_percent", "Mono_int",
                "TNFb_S","CD8")

mtcars <- as.data.frame(estrogen[,which(colnames(estrogen) %in% parameters)])
df <- mtcars[, c("CD8","CXCL9_L")]
df <- na.omit(df)
df$CD8 <- as.numeric(df$CD8)
df$CXCL9_L <- as.numeric(df$CXCL9_L)

a <- ggplot(df, aes(x=CXCL9_L, y=CD8)) + 
  geom_point(color="#791812", size=4)+
  geom_smooth(method="lm")+ylab(expression("T cytotoxic"))+xlab("CXCL9")+
  theme(strip.text.x = element_text(size = 18),
        axis.text.x = element_text(size = 18, angle = 0, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        text = element_text(size = 18))+ggtitle("")+
  theme_bw(base_size = 18)+
  theme(text=element_text(size=18,  family="Arial"))#+expand_limits(y = 0)
c_plot <- a+theme(legend.position = "none")
c_plot
ggsave(filename = paste0("correlation_COPD_CD8_CXCL9.png"), 
       plot=a_plot, width = 4, height = 4, units = c("in"), dpi = 300, scale = 1)

#4)
parameters <- c("DLCOcSB_percent","CXCL9_S", 
                "CXCL10_L","FEV1_FVC_percent",
                "CXCL9_L","CD8", 
                "pCO2_mmHg", "CD8",
                "FEV1_percent", "Mono_int",
                "TNFb_S","CD8")

mtcars <- as.data.frame(estrogen[,which(colnames(estrogen) %in% parameters)])
df <- mtcars[, c("pCO2_mmHg","CD8")]
df <- na.omit(df)
df$CD8 <- as.numeric(df$CD8)
df$pCO2_mmHg <- as.numeric(df$pCO2_mmHg)

a <- ggplot(df, aes(x=CD8, y=pCO2_mmHg)) + 
  geom_point(color="#791812", size=4)+
  geom_smooth(method="lm")+xlab(expression("T cytotoxic"))+ylab(expression("pCO"[2]))+
  theme(strip.text.x = element_text(size = 18),
        axis.text.x = element_text(size = 18, angle = 0, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        text = element_text(size = 18))+ggtitle("")+
  theme_bw(base_size = 18)+
  theme(text=element_text(size=18,  family="Arial"))#+expand_limits(y = 0)
d_plot <- a+theme(legend.position = "none")
d_plot
ggsave(filename = paste0("correlation_COPD_CD8_pCO2.png"), 
       plot=a_plot, width = 4, height = 4, units = c("in"), dpi = 300, scale = 1)


#5)
parameters <- c("DLCOcSB_percent","CXCL9_S", 
                "CXCL10_L","FEV1_FVC_percent",
                "CXCL9_L","CD8", 
                "pCO2_mmHg", "CD8",
                "FEV1_percent", "Mono_int",
                "TNFb_S","CD8")

mtcars <- as.data.frame(estrogen[,which(colnames(estrogen) %in% parameters)])
df <- mtcars[, c("FEV1_percent","Mono_int")]
df <- na.omit(df)
df$FEV1_percent <- as.numeric(df$FEV1_percent)
df$Mono_int <- as.numeric(df$Mono_int)

a <- ggplot(df, aes(x=Mono_int, y=FEV1_percent)) + 
  geom_point(color="#791812", size=4)+
  geom_smooth(method="lm")+xlab(expression("Monocytes intermediate"))+ylab(expression("FEV1%"))+
  theme(strip.text.x = element_text(size = 18),
        axis.text.x = element_text(size = 18, angle = 0, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        text = element_text(size = 18))+ggtitle("")+
  theme_bw(base_size = 18)+
  theme(text=element_text(size=18,  family="Arial"))#+expand_limits(y = 0)
e_plot <- a+theme(legend.position = "none")
e_plot
ggsave(filename = paste0("correlation_COPD_mono_int_FEV.png"), 
       plot=a_plot, width = 4, height = 4, units = c("in"), dpi = 300, scale = 1)

#6)
parameters <- c("DLCOcSB_percent","CXCL9_S", 
                "CXCL10_L","FEV1_FVC_percent",
                "CXCL9_L","CD8", 
                "pCO2_mmHg", "CD8",
                "FEV1_percent", "Mono_int",
                "TNFb_S","CD8")

mtcars <- as.data.frame(estrogen[,which(colnames(estrogen) %in% parameters)])
df <- mtcars[, c("CD8","TNFb_S")]
df <- na.omit(df)
df$CD8 <- as.numeric(df$CD8)
df$TNFb_S <- as.numeric(df$TNFb_S)

a <- ggplot(df, aes(x=TNFb_S, y=CD8)) + 
  geom_point(color="#791812", size=4)+
  geom_smooth(method="lm")+xlab(expression("sTNF-"*beta))+ylab(expression("CD8"^"+"))+
  theme(strip.text.x = element_text(size = 18),
        axis.text.x = element_text(size = 18, angle = 0, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        text = element_text(size = 18))+ggtitle("")+
  theme_bw(base_size = 18)+
  theme(text=element_text(size=18,  family="Arial"))#+expand_limits(y = 0)
f_plot <- a+theme(legend.position = "none")
f_plot
ggsave(filename = paste0("correlation_TNFb_CD8_FEV.png"), 
       plot=a_plot, width = 4, height = 4, units = c("in"), dpi = 300, scale = 1)

c_plot
p2 <- ( b_plot|d_plot|e_plot|c_plot) + patchwork::plot_layout(widths = c(3,3,3,3));p2
ggsave(file = "multipanel_corrplots.png",plot=p2,device = "png",scale=0.9,
       width=20, height=5.2,dpi = 600,units = c("in"))


#save the table ----------
library(openxlsx)
wa <- createWorkbook()
addWorksheet(wa, "network_pval")
addWorksheet(wa, "network_Rcorr")


# Write data to the worksheet (example: write a data frame)
writeData(wa, sheet = "network_pval", dp_sig)
writeData(wa, sheet = "network_Rcorr", x)

saveWorkbook(wa, file = paste0("/home/isilon/users/o_syarif/COPD machine learning/COPD_paper/Fig5/revision_iScience/output/networkdata.xlsx"), overwrite = TRUE)

