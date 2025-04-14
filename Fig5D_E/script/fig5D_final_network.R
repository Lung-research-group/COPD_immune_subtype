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

#directories 
plotdir <- "/home/isilon/users/o_syarif/COPD machine learning/COPD_paper/Fig4H_I/output/plots/"
datainputdir <- "/home/isilon/users/o_syarif/COPD machine learning/COPD_paper/Fig4H_I/input/"
dataoutputdir <- "/home/isilon/users/o_syarif/COPD machine learning/COPD_paper/Fig4H_I/output/data/"

#prepare the master data 
final_data <- read.csv(file =paste0(datainputdir,"final_data_revSC_AvsB.csv"))
elisa_new <- as.data.frame(final_data) # convert to data frame
elisa_new <- elisa_new[which(elisa_new$Consensus_subclass %in% c("A", "B")),] # filtered out everything with no information about subclass
temp <- colnames(elisa_new)# we want to convert to numeric all columns in the data
for (i in temp[!temp %in% c("Diagnosis", "Consensus_subclass")]){
  elisa_new[,i] <- as.numeric(elisa_new[,i])
}

elisa_new$CCL5_S <- NULL #we noticed two parameters are empty (all NAs) and it would generate error in the downstream if we keep them. so we remove them.
elisa_new$TSLP_S <- NULL
elisa_new$NEFA_S <- NULL

clinical <- c("Age","BMI","Smoking_py","CRP",
              "FEV1_percent","FVC_percent","FEV1_FVC_percent",          
              "RV_percent","mPAP","DLCOcSB_percent","pO2_mmHg",                  
              "pCO2_mmHg", "Airspace_Enlargement")
colnames(elisa_new)
sig_parameters <- c("Macs_CD14med_CD1aposHLApos","Macs_CD14hi_CD1aposHLApos",
                    "DC_CD209pos_CD11cposCD1a","DC_CD209neg_CD11cposCd1a",
                    "CD8" , "Mast_CD203pos",
                    "IL1b_S",
                    "GM_CSF_S", "IFNb_S","TNFa_S",
                    "IL10_S",clinical)

up_parameter <- c("Macs_CD14med_CD1aposHLApos","Macs_CD14hi_CD1aposHLApos",
                  "DC_CD209pos_CD11cposCD1a","DC_CD209neg_CD11cposCd1a",
                  "CD8" , "Mast_CD203pos",
                  "IL1b_S",
                  "GM_CSF_S", "IFNb_S","TNFa_S",
                  "IL10_S")
elisa_new <- elisa_new[which(colnames(elisa_new) %in% sig_parameters)]

#prepare the parameter excel file. 
parameter <- read_xlsx(file.choose()) # it should be easily found inside the data folder. Here inside the parameter and the corresponding labels

#calculate fold change (COPD vs DOnor) for each parameter
clinical <- c("Age","BMI","Smoking_py","CRP",
              "FEV1_percent","FVC_percent","FEV1_FVC_percent",          
              "RV_percent","mPAP","DLCOcSB_percent","pO2_mmHg",                  
              "pCO2_mmHg", "Airspace_Enlargement")
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
#write.csv(dat,file = "_COPDvsDonor_revSC_20240216.csv") # save the correlation table, the column "num_of_connection_Var1" later we will remove as it is not representing the true number of connections.
#in the paper, this data frame is available in the worksheet "networks_COPDvsDonor" 

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
#write.csv(fr.all.df,file = paste0("_COPDvsDonorinteractions_revSC.csv")) # we save metadata, in the paper this metadata label in the worksheet "correlations_COPDvsDonor"

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
ggsave(filename = paste0( "rev1_copdsub_sigparameters.png"), 
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
  scale_fill_manual(values=c("darkgreen", "salmon", "orange","lightblue","orchid"), 
                    breaks=c("1","2","3","4","5"))+
  geom_text_repel(data=fr.all.df,aes(x=V1,y=V2,label=expression),parse=T);plot1



#visualise the corr plot--------
corr_var_to <- as.character(dp_sig$to[which(dp_sig$pval <= 0.05)])
corr_var_from <- as.character(dp_sig$from[which(dp_sig$pval <= 0.05)])
main_effects <- unique(c(corr_var_to, corr_var_from))

final_data <- read.csv(file =paste0(datainputdir,"final_data_revSC_AvsB.csv"))
elisa_new <- as.data.frame(final_data) # convert to data frame
elisa_new <- elisa_new[which(elisa_new$Consensus_subclass %in% c("A", "B")),] # filtered out everything with no information about subclass
temp <- colnames(elisa_new)# we want to convert to numeric all columns in the data
for (i in temp[!temp %in% c("Diagnosis", "Consensus_subclass")]){
  elisa_new[,i] <- as.numeric(elisa_new[,i])
}

elisa_new$CCL5_S <- NULL #we noticed two parameters are empty (all NAs) and it would generate error in the downstream if we keep them. so we remove them.
elisa_new$TSLP_S <- NULL
elisa_new$NEFA_S <- NULL

estrogen <- as.matrix(elisa_new[ ,c(main_effects,  "Consensus_subclass")])

#1)
parameters <- c("pO2_mmHg","DC_CD209pos_CD11cposCD1a", 
                "pO2_mmHg","Macs_CD14med_CD1aposHLApos",
                "pO2_mmHg","Macs_CD14hi_CD1aposHLApos", 
                "IL1b_S", "Macs_CD14med_CD1aposHLApos",
                "IL1b_S", "DC_CD209neg_CD11cposCd1a",
                "FEV1_percent","IL1b_S", "Consensus_subclass")

mtcars <- as.data.frame(estrogen[,which(colnames(estrogen) %in% parameters)])
df <- mtcars[, c("pO2_mmHg","DC_CD209pos_CD11cposCD1a","Consensus_subclass")]
df <- na.omit(df)
df$pO2_mmHg <- as.numeric(df$pO2_mmHg)
df$DC_CD209pos_CD11cposCD1a <- as.numeric(df$DC_CD209pos_CD11cposCD1a)

color <- c("#791812",  "#EE4B2B")
names(color) <-  c("A","B")

a <- ggplot(df, aes(x=DC_CD209pos_CD11cposCD1a, y=pO2_mmHg)) + 
  geom_point(aes(color=Consensus_subclass), size=4)+
  geom_smooth(method="lm")+
  ylab(expression("pO"[2]))+xlab(expression("DC CD209"^"+"~"CD11c"^"+"~"CD1a"))+
  theme(strip.text.x = element_text(size = 18),
        axis.text.x = element_text(size = 18, angle = 0, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        text = element_text(size = 18))+ggtitle("")+
  theme_bw(base_size = 18)+
  scale_color_manual(values = color, breaks = names(color))+
  theme(text=element_text(size=18,  family="Arial"))#+expand_limits(y = 0)
a_plot <- a+theme(legend.position = "none")
a_plot

#2)
parameters <- c("pO2_mmHg","DC_CD209pos_CD11cposCD1a", 
                "pO2_mmHg","Macs_CD14med_CD1aposHLApos",
                "pO2_mmHg","Macs_CD14hi_CD1aposHLApos", 
                "IL1b_S", "Macs_CD14med_CD1aposHLApos",
                "IL1b_S", "DC_CD209neg_CD11cposCd1a",
                "FEV1_percent","IL1b_S","Consensus_subclass")

mtcars <- as.data.frame(estrogen[,which(colnames(estrogen) %in% parameters)])
df <- mtcars[, c("pO2_mmHg","Macs_CD14med_CD1aposHLApos","Consensus_subclass")]
df <- na.omit(df)
df$pO2_mmHg <- as.numeric(df$pO2_mmHg)
df$Macs_CD14med_CD1aposHLApos <- as.numeric(df$Macs_CD14med_CD1aposHLApos)

a <- ggplot(df, aes(x=Macs_CD14med_CD1aposHLApos, y=pO2_mmHg)) + 
  geom_point(aes(color=Consensus_subclass), size=4)+
  geom_smooth(method="lm")+ylab(expression("pO"[2]))+xlab(expression("Macs CD14"^"med"~"CD1a"^"+"~"HLA"^"+"))+
  theme(strip.text.x = element_text(size = 18),
        axis.text.x = element_text(size = 18, angle = 0, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        text = element_text(size = 18))+ggtitle("")+
  theme_bw(base_size = 18)+
  scale_color_manual(values = color, breaks = names(color))+
  theme(text=element_text(size=18,  family="Arial"))#+expand_limits(y = 0)
b_plot <- a+theme(legend.position = "none")
b_plot


#3)
parameters <- c("pO2_mmHg","DC_CD209pos_CD11cposCD1a", 
                "pO2_mmHg","Macs_CD14med_CD1aposHLApos",
                "pO2_mmHg","Macs_CD14hi_CD1aposHLApos", 
                "IL1b_S", "Macs_CD14med_CD1aposHLApos",
                "IL1b_S", "DC_CD209neg_CD11cposCd1a",
                "FEV1_percent","IL1b_S","Consensus_subclass")

mtcars <- as.data.frame(estrogen[,which(colnames(estrogen) %in% parameters)])
df <- mtcars[, c("pO2_mmHg","Macs_CD14hi_CD1aposHLApos","Consensus_subclass")]
df <- na.omit(df)
df$Macs_CD14hi_CD1aposHLApos <- as.numeric(df$Macs_CD14hi_CD1aposHLApos)
df$pO2_mmHg <- as.numeric(df$pO2_mmHg)

a <- ggplot(df, aes(x=Macs_CD14hi_CD1aposHLApos, y=pO2_mmHg)) + 
  geom_point(aes(color=Consensus_subclass), size=4)+
  geom_smooth(method="lm")+ylab(expression("pO"[2]))+xlab("Macs CD14"^"hi"~"CD1a"^"+"~"HLA"^"+")+
  theme(strip.text.x = element_text(size = 18),
        axis.text.x = element_text(size = 18, angle = 0, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        text = element_text(size = 18))+ggtitle("")+
  theme_bw(base_size = 18)+
  scale_color_manual(values = color, breaks = names(color))+
  theme(text=element_text(size=18,  family="Arial"))#+expand_limits(y = 0)
c_plot <- a+theme(legend.position = "none")
c_plot

#4)
parameters <- c("pO2_mmHg","DC_CD209pos_CD11cposCD1a", 
                "pO2_mmHg","Macs_CD14med_CD1aposHLApos",
                "pO2_mmHg","Macs_CD14hi_CD1aposHLApos", 
                "IL1b_S", "Macs_CD14med_CD1aposHLApos",
                "IL1b_S", "DC_CD209neg_CD11cposCd1a",
                "FEV1_percent","IL1b_S","Consensus_subclass")

mtcars <- as.data.frame(estrogen[,which(colnames(estrogen) %in% parameters)])
df <- mtcars[, c("IL1b_S","Macs_CD14med_CD1aposHLApos","Consensus_subclass")]
df <- na.omit(df)
df$Macs_CD14med_CD1aposHLApos <- as.numeric(df$Macs_CD14med_CD1aposHLApos)
df$IL1b_S <- as.numeric(df$IL1b_S)

a <- ggplot(df, aes(x=Macs_CD14med_CD1aposHLApos, y=IL1b_S)) + 
  geom_point(aes(color=Consensus_subclass), size=4)+
  geom_smooth(method="lm")+xlab(expression("Macs CD14"^"med"~"CD1a"^"+"~"HLA"^"+"))+ylab(expression("sIL-1"*beta))+
  theme(strip.text.x = element_text(size = 18),
        axis.text.x = element_text(size = 18, angle = 0, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        text = element_text(size = 18))+ggtitle("")+
  theme_bw(base_size = 18)+
  scale_color_manual(values = color, breaks = names(color))+
  theme(text=element_text(size=18,  family="Arial"))#+expand_limits(y = 0)
d_plot <- a+theme(legend.position = "none")
d_plot


#5)
parameters <- c("pO2_mmHg","DC_CD209pos_CD11cposCD1a", 
                "pO2_mmHg","Macs_CD14med_CD1aposHLApos",
                "pO2_mmHg","Macs_CD14hi_CD1aposHLApos", 
                "IL1b_S", "Macs_CD14med_CD1aposHLApos",
                "IL1b_S", "DC_CD209neg_CD11cposCd1a",
                "FEV1_percent","IL1b_S","Consensus_subclass")

mtcars <- as.data.frame(estrogen[,which(colnames(estrogen) %in% parameters)])
df <- mtcars[, c("IL1b_S","DC_CD209neg_CD11cposCd1a","Consensus_subclass")]
df <- na.omit(df)
df$IL1b_S <- as.numeric(df$IL1b_S)
df$DC_CD209neg_CD11cposCd1a <- as.numeric(df$DC_CD209neg_CD11cposCd1a)

a <- ggplot(df, aes(x=DC_CD209neg_CD11cposCd1a, y=IL1b_S)) + 
  geom_point(aes(color=Consensus_subclass), size=4)+
  geom_smooth(method="lm")+xlab(expression("DC CD209"^"-"~"CD11c"^"+"~"CD1a"))+ylab(expression("sIL-1"*beta))+
  theme(strip.text.x = element_text(size = 18),
        axis.text.x = element_text(size = 18, angle = 0, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        text = element_text(size = 18))+ggtitle("")+
  theme_bw(base_size = 18)+
  scale_color_manual(values = color, breaks = names(color))+
  theme(text=element_text(size=18,  family="Arial"))#+expand_limits(y = 0)
e_plot <- a+theme(legend.position = "none")
e_plot


#6)
parameters <- c("pO2_mmHg","DC_CD209pos_CD11cposCD1a", 
                "pO2_mmHg","Macs_CD14med_CD1aposHLApos",
                "pO2_mmHg","Macs_CD14hi_CD1aposHLApos", 
                "IL1b_S", "Macs_CD14med_CD1aposHLApos",
                "IL1b_S", "DC_CD209neg_CD11cposCd1a",
                "FEV1_percent","IL1b_S","Consensus_subclass")

mtcars <- as.data.frame(estrogen[,which(colnames(estrogen) %in% parameters)])
df <- mtcars[, c("FEV1_percent","IL1b_S","Consensus_subclass")]
df <- na.omit(df)
df$IL1b_S <- as.numeric(df$IL1b_S)
df$FEV1_percent <- as.numeric(df$FEV1_percent)

a <- ggplot(df, aes(x=IL1b_S, y=FEV1_percent)) + 
  geom_point(aes(color=Consensus_subclass), size=4)+
  geom_smooth(method="lm")+xlab(expression("sIL-1"*beta))+ylab(expression("FEV1%"))+
  theme(strip.text.x = element_text(size = 18),
        axis.text.x = element_text(size = 18, angle = 0, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        text = element_text(size = 18))+ggtitle("")+
  theme_bw(base_size = 18)+
  scale_color_manual(values = color, breaks = names(color))+
  theme(text=element_text(size=18,  family="Arial"))#+expand_limits(y = 0)
f_plot <- a+theme(legend.position = "none")
f_plot


c_plot
p2 <- ( a_plot|b_plot|c_plot|d_plot|e_plot|f_plot) + 
  patchwork::plot_layout(widths = c(3,3,3,3,3,3));p2
ggsave(file = "multipanel_COPDsub.png",plot=p2,device = "png",scale=0.9,
       width=30, height=5.2,dpi = 600,units = c("in"))


#save the table ----------
library(openxlsx)
wa <- createWorkbook()
addWorksheet(wa, "network_pval")
addWorksheet(wa, "network_Rcorr")


# Write data to the worksheet (example: write a data frame)
writeData(wa, sheet = "network_pval", dp_sig)
writeData(wa, sheet = "network_Rcorr", x)

saveWorkbook(wa, file = paste0("/home/isilon/users/o_syarif/COPD machine learning/COPD_paper/Fig5/revision_iScience/output/networkdata_subtype.xlsx"), overwrite = TRUE)

