#will perform GMM and kmeans clustering with set to 2 clusters, exports results
#---Load libraries:########--------#####
	library(readxl)     ##read in excel file 
  library(missMDA)    ##imputePCA to impute missings in dataset to enable later classical PCA
  library(ggplot2)    ##for plotting 
	library(ggpubr)     ##for plotting 
  library(ggrepel)    ##for plotting labels not on top of each other
  library(colorspace) ##darken/lighten
	library(factoextra) ##for PCA loadings plots with arrows 
  library(mclust)     ##for GMM and kmeans clustering

  theme_set(theme_pubr(base_size=6, border = T, legend = "none"))
  theme_journal <- theme(panel.grid.major = element_line(color = "light grey", linewidth = 0.2, linetype = "solid"), 
                         panel.grid.minor = element_line(color = "light grey", linewidth = 0.05, linetype = "solid"), 
                         plot.title = element_text(size = rel(1)),
                         axis.line.x = element_line(color="black", linewidth = 0.2),
                         axis.line.y = element_line(color="black", linewidth = 0.2),
                         axis.text = element_text(color="black", size = 6),
                         panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.2)
  )
  
  w = 49; h = 52;  #size settings for plot exports

#---Initializing:  setting all names, sizes, colors, read in data etc.: #####			
	#Set path where to find all .R files and folders with input/output
	statisticsfilepath <- paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/")
	model_name_per <- "15_hFACS_COPD_Lung_C_per_LOG_x1c_filtDD"
	model_name_cc <- "16_hFACS_COPD_Lung_C_cc_LOG_x1c_filtDD"
	input_file_per <- paste0(statisticsfilepath,"input data/",model_name_per,".xlsx")
	input_file_cc <- paste0(statisticsfilepath,"input data/",model_name_cc,".xlsx")
	output_folder <- paste0(statisticsfilepath,"output_Fig4_COPD_subclass/")
	
	#read in data	
	unique_sample_ID_colname <- "Unique_Sample_ID" #column with unique identifier per row in that dataset
	subdataset_per <- data.frame(read_excel(input_file_per,  sheet = "Sheet 1"))
	sample_info_per <- subdataset_per[,1:3]
	rownames(subdataset_per) <- rownames(sample_info_per) <- subdataset_per[,unique_sample_ID_colname]
	subdataset_per <- subdataset_per[,-c(1:3)]

	subdataset_cc <- data.frame(read_excel(input_file_cc,  sheet = "Sheet 1"))
	sample_info_cc <- subdataset_cc[,1:3]
	rownames(subdataset_cc) <- rownames(sample_info_cc) <- subdataset_cc[,unique_sample_ID_colname]
	subdataset_cc <- subdataset_cc[,-c(1:3)]

#---Impute missings and perform classical PCA: ######
	if(sum(is.na(subdataset_per))!=0) { #impute if missings exist
	  if(ncol(subdataset_per)<=10) {ncomp_impute <- ncol(subdataset_per)-1} else {ncomp_impute <- 10}
	  subdataset_imputed_per <- imputePCA(subdataset_per,ncomp_impute)$completeObs   
	} else {
	  subdataset_imputed_per <- subdataset_per
	}
	
	if(sum(is.na(subdataset_cc))!=0) { #impute if missings exist
	  if(ncol(subdataset_cc)<=10) {ncomp_impute <- ncol(subdataset_cc)-1} else {ncomp_impute <- 10}
	  subdataset_imputed_cc <- imputePCA(subdataset_cc,ncomp_impute)$completeObs   
	} else {
	  subdataset_imputed_cc <- subdataset_cc
	}
	
#---Perform GMM and kmeans Clustering and export results: ######			    	
 
	#to allow matching to our naming of cluster results
	subclass <- as.numeric(gsub("B", "2", gsub("A", "1", sample_info_per[,"COPD_subclass"]))) #same as in _cc
  names(subclass) <- sample_info_per[,unique_sample_ID_colname]
  
  set.seed(234971) #so that kmeans has always same result
    
  gmm_res_per <- Mclust(subdataset_imputed_per, G = 2) #G number of clusters
  gmm_plot <- fviz_cluster(gmm_res_per, geom = "point", data = subdataset_imputed_per, 
                           ellipse.type = "convex", main = "GMM Clustering into 2") + 
      geom_text_repel(aes(label = rownames(subdataset_imputed_per))) + 
    labs(caption = paste0(model_name_per,"_gmm_2clust.png")) + theme_journal
  ggsave(paste0(output_folder, model_name_per,"_gmm_2clust.png"), 
         plot = gmm_plot, dpi = 600, width = w*3, height = h*3, units = "mm"); #dev.off()
  Diff_gmm_per <- gmm_res_per$classification == subclass
  
  gmm_res_cc <- Mclust(subdataset_imputed_cc, G = 2) #G number of clusters
  gmm_plot <- fviz_cluster(gmm_res_cc, geom = "point", data = subdataset_imputed_cc, 
                           ellipse.type = "convex", main = "GMM Clustering into 2") + 
    geom_text_repel(aes(label = rownames(subdataset_imputed_per))) + 
    labs(caption = paste0(model_name_cc,"_gmm_2clust.png")) + theme_journal
  ggsave(paste0(output_folder, model_name_cc,"_gmm_2clust.png"), 
         plot = gmm_plot, dpi = 600, width = w*3, height = h*3, units = "mm"); #dev.off()
  Diff_gmm_cc <- gmm_res_cc$classification == subclass
    
  km_res_per <- kmeans(subdataset_imputed_per, centers = 2, nstart = 5)
  km_plot <- fviz_cluster(km_res_per, geom = "point", data = subdataset_imputed_per, 
                          ellipse.type = "convex", main = "K-means Clustering into 2") + 
      geom_text_repel(aes(label = rownames(subdataset_imputed_per))) +
    labs(caption = paste0(model_name_per,"_kmeans_2clust.png")) + theme_journal
  ggsave(paste0(output_folder, model_name_per,"_kmeans_2clust.png"), 
         plot = km_plot, dpi = 600, width = w*3, height = h*3, units = "mm"); #dev.off()
  Diff_km_per <- km_res_per$cluster == subclass
  
  km_res_cc <- kmeans(subdataset_imputed_cc, centers = 2, nstart = 5)
  km_plot <- fviz_cluster(km_res_cc, geom = "point", data = subdataset_imputed_cc, 
                          ellipse.type = "convex", main = "K-means Clustering into 2") + 
    geom_text_repel(aes(label = rownames(subdataset_imputed_cc))) +
    labs(caption = paste0(model_name_cc,"_kmeans_2clust.png")) + theme_journal
  ggsave(paste0(output_folder, model_name_cc,"_kmeans_2clust.png"), 
         plot = km_plot, dpi = 600, width = w*3, height = h*3, units = "mm"); #dev.off()
  Diff_km_cc <- km_res_cc$cluster == subclass  
    
 clusterids <- data.frame(Unique_Sample_ID = names(subclass), 
                          km_res_per = km_res_per$cluster, 
                          gmm_res_per = gmm_res_per$classification, 
                          km_res_cc = km_res_cc$cluster, 
                          gmm_res_cc = gmm_res_cc$classification, 
                          COPD_subclass = sample_info_per[,"COPD_subclass"],
                          COPD_subclass_Nr = subclass, 
                          Diff_km_per = Diff_km_per, 
                          Diff_gmm_per = Diff_gmm_per,
                          Diff_km_cc = Diff_km_cc, 
                          Diff_gmm_cc = Diff_gmm_cc)
 
 write.table(clusterids, file = paste0(output_folder,"Fig4_COPD_subclass_GMM_kmeans_clustering.txt"), 
                sep = "\t", row.names = FALSE, col.names = TRUE)
