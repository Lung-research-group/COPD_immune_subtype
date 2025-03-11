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

#calculate the significant parameters
storage <- data.frame() #we want to calculate significance in every parameters using wilcoxon test. There should not be any NAs. 
for(i in names(elisa_new)[!names(elisa_new) %in% c("Consensus_subclass", "Diagnosis", "Unique_Sample_ID")]){ # we exclude column "Consensus_subclass" and "Unique Sample ID" from the calculations. 
  elisa <- elisa_new[which(!is.na(elisa_new[, i])),]  # we remove the NAs
  x <- compare_means(eval(parse(text=paste0(i, " ~ Consensus_subclass"))), data=elisa, p.adjust.method = "fdr") #using compare means function from ggpubr to calculate the wilcoxon test (A vs B)
  storage <- rbind(storage, x) # then bind all results to one column
}

#prepare the parameter excel file. 
parameter <- read_xlsx(path=paste0(datainputdir, "parameter_rev01.xlsx")) # it should be easily found inside the data folder. Here inside the parameter and the corresponding labels

#calculate the effect size
x <-lapply(elisa_new[!names(elisa_new)  %in% c("Consensus_subclass", "Diagnosis", "Unique_Sample_ID")],function(x) cohen.d(x ~ elisa_new$Consensus_subclass))# we exclude column "Consensus_subclass" and "Unique Sample ID" from the calculations.  inside we can find list of cohen's caculations for all parameters
cohen_d <- c()  # create new array that we will use to store all results.
for(i in names(elisa_new)[!names(elisa_new) %in% c("Consensus_subclass", "Diagnosis", "Unique_Sample_ID")]){ #we always exclude column "Consensus_subclass" from the calculations. 
  a <- x[[i]][["estimate"]] # we take estimate from each parameter's calculation
  cohen_d[i] <- a #and call it as "a" elemnt in cohen_d array
}
temp <- as.data.frame(cohen_d) # we convert to data frame
temp$.y. <- rownames(temp) #create a new column

#merge significance and effect size calculations
storage <- merge(storage, temp, by=".y.") #we merge the cohen's values with significance values from previous calculations

#calculate fold change (A vs B) for each parameter
clinical <- c("Age","BMI","Smoking_py","CRP", "Airspace_Enlargement",
              "FEV1_percent","FVC_percent","FEV1_FVC_percent",          
              "RV_percent","mPAP","DLCOcSB_percent","pO2_mmHg",                  
              "pCO2_mmHg")
cells <- c( "Basophils", "CD19","CD3",                       
            "CD4","CD8","DC_CD209neg_CD11cposCd1a",
            "DC_CD209pos_CD11cposCD1a","gdTCR","Macs_CD14hi_CD1anegHLApos",
            "Macs_CD14hi_CD1aposHLAneg",  "Macs_CD14hi_CD1aposHLApos", 
            "Macs_CD14med_CD1anegHLApos","Macs_CD14med_CD1aposHLAneg", 
            "Macs_CD14med_CD1aposHLApos","Monocytes","Neutrophils", "Mast_CD203neg","Mast_CD203pos",
            "NKcells", "NKT","pDC","Tregs") 
cytokine <- colnames(final_data)[-which(colnames(final_data) %in% c(clinical, cells, "Unique_Sample_ID","Diagnosis", "Consensus_subclass"))]

elisa_new$Unique_Sample_ID <- NULL  # we remove the column "Unique Sample ID" from the data
melt <- melt(elisa_new)# melt the data so we have column diagnosis, varibale (parameter), and values
x_cells_cytokine <-melt[which(melt$variable %in% c(cells, cytokine)),] %>% na.omit()%>% # we group by variable and calculate the fold change between A and B
  group_by(variable)  %>%
  summarise(fold_change = mean(value[Consensus_subclass=="A"], na.rm = T)- mean(value[Consensus_subclass=="B"], na.rm = T))

x_clinical <-melt[which(melt$variable %in% c(clinical)),] %>% na.omit()%>% # we group by variable and calculate the fold change between A and B
  group_by(variable)  %>%
  summarise(fold_change = mean(value[Consensus_subclass=="A"], na.rm = T)/ mean(value[Consensus_subclass=="B"], na.rm = T))

x_cells_cytokine$direction <- case_when(x_cells_cytokine$fold_change >= 0 ~ "upregulated",
                                        x_cells_cytokine$fold_change < 0 ~ "downregulated") # annotate the direction based on the fold change information

x_clinical$direction <- case_when(x_clinical$fold_change >= 1 ~ "upregulated",
                                  x_clinical$fold_change < 1 ~ "downregulated") # annotate the direction based on the fold change information

x <- rbind(x_clinical, x_cells_cytokine)
#merge the significance, effect size, and fold change data
storage <- merge(storage, x, by.x=".y.", by.y="variable")

#visualisation
#create the igraph
main.effects <- storage$.y.[which(storage$p.adj <0.1)] # we only select significant parameters to be visualised
estrogen <- as.matrix(elisa_new[ ,main.effects])# we filter the data only  keeping the significant parameters

set.seed(seed = 123)
g <- graph.adjacency( # Create a graph adjacency based on correlation distances between parameters in  pairwise fashion.
  as.matrix(as.dist(cor(estrogen, method="spearman", use="pairwise.complete.obs"))), # so basically those which correlated highly will be positioned closer in the graph
  mode="undirected",
  weighted=TRUE,# we want to weigh the edge (weight here is the correlation coefficient)
  diag=FALSE
)

g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE) # Simplfy the adjacency object

g <- delete_edges(g, E(g)[which(abs(E(g)$weight)<0.5)])# Remove edges below absolute Pearson correlation 0.5 and with NA value
g <- delete_edges(g, E(g)[!complete.cases(E(g)$weight)])

V(g)$name <- V(g)$name # Assign names to the graph vertices (optional). vertices here is the parameters name

y <- get.data.frame(g) # get the igraph object with weight of each connection

y$change <- case_when(y$weight > 0 ~ "positive corr", # annotate the weight or the correlation coeficient
                      y$weight < 0 ~ "negative corr")
y$corr <- paste0(y$from, "_", y$to) # label the connecting parameters
E(g)$weight <- abs(E(g)$weight) # add the weight into edges information. Edges here is the connecting parameters

#create metadata
fr.all <- layout.fruchterman.reingold(g) # we apply fruchterman reingold to draw the graph
fr.all.df <- as.data.frame(fr.all) # we get the coordinate x and y as V1 and V2 columns
fr.all.df$species <- colnames(estrogen) #we attach parameters name into coordinate values
fr.all.df <- merge(fr.all.df, storage, by.x="species", by.y=".y.")# we attach the rest of the data (significance values, foldchange, effect size/cohen values) to the coordinate data frame
fr.all.df$abs_cohen <- abs(fr.all.df$cohen_d)# we convert cohen's values to absolute values as we will use to draw the size of the node
fr.all.df <- merge(fr.all.df, parameter, by.x="species", by.y="Parameter_Name")# we attach the parameter label (for official COPD paper annotation)


#create correlation data frame to be saved as supplementary
x <- get.data.frame(g)
x$corr <- paste0(x$from, "_", x$to)
x$change <- y$change[match(x$corr, y$corr)]#attach the information from y table above, the information about positiive and negative corr
x$from.x <- fr.all.df$V1[match(x$from, fr.all.df$species)]*20 #we took V1 parameters from metadata above, and put them here in table as from.x column;match with paramater in  column "from" to  column "species" in metadata for coordinate x (V1)
x$from.y <- fr.all.df$V2[match(x$from, fr.all.df$species)]*20 # match with paramater in  column "from" to  column "species" in metadata for coordinate y (V2)

x$to.x <- fr.all.df$V1[match(x$to, fr.all.df$species)]*20   #match with paramater in  column "to" to  column "species" in metadata for coordinate y (V1)
x$to.y <- fr.all.df$V2[match(x$to, fr.all.df$species)]*20 #match with paramater in  column "to" to  column "species" in metadata for coordinate y (V2)


fr.all.df$V1 <- fr.all.df$V1*20 # copy the coordinate and increase the coordinate as we did previously to other connections
fr.all.df$V2 <- fr.all.df$V2*20


x$change <- relevel(as.factor(x$change), ref="positive corr") # we relevel the factor of correlation label
fr.all.df$label <- ordered(fr.all.df$label, levels=c("clinical and lung function","circulation","lung"))  # we also reorder the level of origin label of parameters

temp <- as.data.frame(table(x$from))  # we calculate connections "from". we have two connections: from and to
colnames(temp)[2] <- "num_of_connection_Var1" # we label the name of connections
dat <- y[which(y$from %in% temp$Var1),] # filter data, only including the confirmed connections from both x and y (just as double check)
dat <- merge(dat, temp, by.x="from", by.y="Var1")# attach the rest information

write.csv(dat,file = paste0(dataoutputdir,"_AvsB_revSC.csv"))# save the correlation table, the column "num_of_connection_Var1" later we will remove as it is not representing the true number of connections.
#in the paper, this data frame is available in the worksheet "networks_AvsB" 

#cluster the parameters using the fast greedy algorithm
fc <- cluster_fast_greedy(g)# we take the igraph object #cluster_fast_greedy #cluster_optimal #cluster_leading_eigen
cluster <- cut_at(fc, 2)
mem <-membership(fc)# assign clusters
z <- names(mem)#take all parameters name
y <- as.numeric(cluster)# take all cluster number
h <-data.frame(name =z, community =y)# make data frame of cluster name and parameters name
fr.all.df <- merge(fr.all.df, h, by.x="species", by.y="name")# attach the cluster name and parameters name to metadata we created above

a <- as.data.frame(table(dat$to)) # we take parameters name in column "to" from correlation data frame, and calculate how many times they appear indicating how many connections they have as from/origin
b <- as.data.frame(table(dat$from)) # we take parameters name in column "from" from correlation data frame, and calculate how many times they appear indicating how many connections they have as to/destination
hh <- full_join(a,b, by="Var1")# we join the information of above a and b
hh[is.na(hh)] <- 0# we assign any NAs to 0
hh$total <- hh$Freq.x+hh$Freq.y # we sum the total number of connections
fr.all.df <- merge(fr.all.df, hh, by.x="species", by.y="Var1") # and attach this information to metadata

write.csv(fr.all.df,file = paste0(dataoutputdir,"_AvsBinteractions_revSC.csv"))  # we save metadata, in the paper this metadata label in the worksheet "correlations_COPDvsDonor"

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
  text = element_text(size=20),
  panel.grid.major = element_blank(),  #remove major-grid labels
  panel.grid.minor = element_blank())  #remove minor-grid labels

plot1 <- ggplot() +
  geom_link(data=x,aes(x=from.x,xend = to.x, y=from.y,yend = to.y, alpha=weight,size=weight/15, linetype=change),color="black") + # add line type
  geom_point(data=fr.all.df,aes(x=V1,y=V2,fill=as.factor(community), size=abs_cohen, shape=direction),show_guide=T)+ # add colour scaling for group membership
  scale_shape_manual(values=c(25,24))+
  theme_set(theme_pubr(base_size=20, border = T))+
  scale_alpha(range = c(0.4, 0.8))+
  scale_size(range = c(0,10))+
  scale_x_continuous(expand=c(0,4))+  # expand the x limits 
  scale_y_continuous(expand=c(0,4))+ # expand the y limits
  theme_bw()+  # use the ggplot black and white theme
  theme_journal+
  geom_text_repel(data=fr.all.df,aes(x=V1,y=V2,label=expression, point.size=abs_cohen*2),parse=T)
plot1

ggsave(filename = paste0(plotdir,"cutoff0.5edge__AvsB01community_revSC_02.png"), plot=plot1, width = 10, height = 8, 
       units = c("in"), dpi = 300, scale = 0.8)

#with parameter labels
plot1 <- ggplot() +
  geom_link(data=x,aes(x=from.x,xend = to.x, y=from.y,yend = to.y, alpha=weight,size=weight/30, linetype=change),color="black") + # add line type
  geom_point(data=fr.all.df,aes(x=V1,y=V2,fill=as.factor(label), size=abs_cohen, shape=direction),show_guide=T)+ # add colour scaling for group membership
  scale_shape_manual(values=c(25,24))+
  theme_set(theme_pubr(base_size=20, border = T))+
  scale_alpha(range = c(0.4, 0.8))+
  scale_fill_manual(values=c("red", "royalblue4", "goldenrod"), 
                    breaks=c("circulation","lung","clinical and lung function"))+
  scale_size(range = c(0,10))+
  scale_x_continuous(expand=c(0,4))+  # expand the x limits 
  scale_y_continuous(expand=c(0,4))+ # expand the y limits
  theme_bw()+  # use the ggplot black and white theme
  theme_journal+
  geom_text_repel(data=fr.all.df,aes(x=V1,y=V2,label=expression, point.size=abs_cohen*2),parse=T)
plot1
ggsave(filename = paste0(plotdir, "cutoff0.5edge__AvsB01_revSC_02.png"), plot=plot1, 
       width = 10, height = 8, units = c("in"), dpi = 300, scale = 0.8)

#legend!
plot1 <- ggplot() +
  geom_link(data=x,aes(x=from.x,xend = to.x, y=from.y,yend = to.y, alpha=weight,size=weight/30, linetype=change),color="black") + # add line type
  geom_point(data=fr.all.df,aes(x=V1,y=V2,fill=as.factor(label), size=abs_cohen, shape=direction),show_guide=T)+ # add colour scaling for group membership
  scale_shape_manual(values=c(25,24))+
  theme_set(theme_pubr(base_size=20, border = T))+
  scale_alpha(range = c(0.4, 0.8))+
  scale_fill_manual(values=c("red", "royalblue4", "goldenrod"), 
                    breaks=c("circulation","lung","clinical and lung function"))+
  scale_size(range = c(0,10))+
  scale_x_continuous(expand=c(0,4))+  # expand the x limits 
  scale_y_continuous(expand=c(0,4))+ # expand the y limits
  theme_bw()+  # use the ggplot black and white theme
  theme_journal+theme(legend.position = "right")+
  geom_text_repel(data=fr.all.df,aes(x=V1,y=V2,label=expression, point.size=abs_cohen*2),parse=T)
plot1
ggsave(filename = paste0(plotdir, "legend_cutoff0.5edge__COPDvsDonor01_revSC.png"), 
       plot=plot1, width = 10, height = 10, units = c("in"), dpi = 300, scale = 0.8)

ggplot() + # add line type
  geom_point(data=fr.all.df,aes(x=V1,y=V2,fill=as.factor(label), size=abs_cohen, shape=direction),show_guide=T)+ # add colour scaling for group membership
  scale_shape_manual(values=c(25,24))+
  theme_set(theme_pubr(base_size=20, border = T))+
  scale_alpha(range = c(0.4, 0.8))+
  scale_fill_manual(values=c("red", "royalblue4", "goldenrod"), 
                    breaks=c("circulation","lung","clinical and lung function"))+
  scale_size(range = c(0,10))+
  scale_x_continuous(expand=c(0,4))+  # expand the x limits 
  scale_y_continuous(expand=c(0,4))+ # expand the y limits
  theme_bw()+  # use the ggplot black and white theme
  theme_journal+theme(legend.position = "right")+
  geom_text_repel(data=fr.all.df,aes(x=V1,y=V2,label=expression, point.size=abs_cohen*2),parse=T)
ggsave(filename = paste0(plotdir, "cutoff0.5edge__AvsB01_revSC_legend_point.png"),
       width = 12, height = 8, units = c("in"), dpi = 300, scale = 0.8)


#correlation plots for selected parameters
#visualisation
#create the igraph
main.effects <- storage$.y.[which(storage$p.adj <0.1)] # we only select significant parameters to be visualised
estrogen <- as.matrix(elisa_new[ ,main.effects]) # we filter the data only  keeping the significant parameters

cor.mtest <- function(mat) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      error <- try(tmp <- cor.test(mat[, i], mat[, j]),
                   silent =T)
      if (class(error) == "try-error") {
        p.mat[i, j] <- NA
      } else {
        p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      }
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

c <- cor.mtest(estrogen)
cor <- as.matrix(cor(estrogen, method="spearman", use="pairwise.complete.obs"))
x <- c

c[] <- p.adjust(c, method = "fdr")

melt <- melt(cor)
melt_pval <- melt(c)
melt$variable <- paste0(melt$Var1, "_", melt$Var2)
colnames(melt)[3] <- "fold_change"
colnames(melt_pval)[3] <- "adjusted_pvalue"
melt_pval$variable <- paste0(melt_pval$Var1, "_", melt_pval$Var2)


table_correlation <-merge(melt, melt_pval, by="variable")
table_correlation <- table_correlation[, c("fold_change", "adjusted_pvalue", "variable")]
colnames(table_correlation)[1] <- "correlation"
table_correlation$correlation <- as.numeric(table_correlation$correlation)
corrrelation <- table_correlation[which(abs(table_correlation$adjusted_pvalue) >0),]
corrrelation$abs <- abs(corrrelation$correlation)

corrrelation <- corrrelation[which(corrrelation$variable %in% dat$corr),]

#visualise the corr plots
main.effects <- storage$.y.[which(storage$p.adj <0.1)] # we only select significant parameters to be visualised
estrogen <- as.matrix(elisa_new[ ,c(main.effects,"Consensus_subclass")]) # we filter the data only  keeping the significant parameters
parameters <- c("DC_CD209pos_CD11cposCD1a","pO2_mmHg",
                "TNFa_S",
                "Macs_CD14med_CD1aposHLApos", "GM_CSF_S", "Consensus_subclass")

mtcars <- as.data.frame(estrogen[,which(colnames(estrogen) %in% parameters)])

df <- mtcars[, c("DC_CD209pos_CD11cposCD1a","pO2_mmHg", "Consensus_subclass")]
df <- na.omit(df)
df$DC_CD209pos_CD11cposCD1a <- as.numeric(df$DC_CD209pos_CD11cposCD1a)
df$pO2_mmHg <- as.numeric(df$pO2_mmHg)

color <- c("#791812",  "#EE4B2B")
names(color) <-  c("A","B")

a <- ggplot(df, aes(x=DC_CD209pos_CD11cposCD1a, y=pO2_mmHg)) + 
  geom_point(aes(color=Consensus_subclass), size=4)+
  geom_smooth(method="lm")+
  theme(strip.text.x = element_text(size = 16),
        axis.text.x = element_text(size = 16, angle = 0, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        text = element_text(size = 16))+ggtitle("")+
  theme_bw(base_size = 16)+xlab("DC CD209"^"+"~"CD11c"^"+"~"CD1a")+ylab("pO2")+
  scale_color_manual(values = color, breaks = names(color))+
  theme(text=element_text(size=16,  family="Arial"))
a
a_plot <- a+theme(legend.position = "none")
a_plot

df <- mtcars[, c("DC_CD209pos_CD11cposCD1a","TNFa_S", "Consensus_subclass")]
df <- na.omit(df)
df$DC_CD209pos_CD11cposCD1a<- as.numeric(df$DC_CD209pos_CD11cposCD1a)
df$TNFa_S <- as.numeric(df$TNFa_S)
b <- ggplot(df, aes(x=DC_CD209pos_CD11cposCD1a, y=TNFa_S)) + 
  geom_point(aes(color=Consensus_subclass), size=4)+
  geom_smooth(method="lm")+
  theme(strip.text.x = element_text(size = 16),
        axis.text.x = element_text(size = 16, angle = 0, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        text = element_text(size = 16))+ggtitle("")+
  theme_bw(base_size = 16)+xlab("DC CD209"^"+"~"CD11c"^"+"~"CD1a")+ylab("TNF-α")+
  scale_color_manual(values = color, breaks = names(color))+
  theme(text=element_text(size=16,  family="Arial"))
b
b_plot <- b+theme(legend.position = "none")
b_plot


df <- mtcars[, c("Macs_CD14med_CD1aposHLApos","GM_CSF_S", "Consensus_subclass")]
df <- na.omit(df)
df$Macs_CD14med_CD1aposHLApos <- as.numeric(df$Macs_CD14med_CD1aposHLApos)
df$GM_CSF_S <- as.numeric(df$GM_CSF_S)
c <- ggplot(df, aes(x=Macs_CD14med_CD1aposHLApos, y=GM_CSF_S)) + 
  geom_point(aes(color=Consensus_subclass), size=4)+
  geom_smooth(method="lm")+
  theme(strip.text.x = element_text(size = 16),
        axis.text.x = element_text(size = 16, angle = 0, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        text = element_text(size = 16))+ggtitle("")+
  theme_bw(base_size = 16)+ylab("sGM-CSF")+xlab("Macs CD14"^"med"~"CD1a"^"+"~"HLA"^"+")+
  scale_color_manual(values = color, breaks = names(color))+
  theme(text=element_text(size=16,  family="Arial"))
c
c_plot <- c+theme(legend.position = "none")
c_plot

ggsave(filename = paste0(plotdir, "correlation_AvsB_b.png"), 
       plot=b_plot, width = 5, height = 5, units = c("in"), dpi = 300, scale = 1)
correlation_AvsB <- corrrelation
correlation_AvsB$abs <- abs(correlation_AvsB$correlation)

write.csv(correlation_AvsB, file =paste0(plotdir, "_correlation_AvsB_table.csv"))


