#Creates a Heatmap of the data (z-scaled)
#---Load libraries:########--------#####
	library(readxl)       #read in excel file 
  library(ggplot2)      #for plotting 
  library(pheatmap)     #heatmap itsels
  library(RColorBrewer) #make color scale gradients
  library(dendsort)     #sort the dendograms, increases visibility of clustering results
  library(dplyr) #left_join()

  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))  #sorting helper function

#---Initializing:  setting all names, sizes, colors, read in data etc.: #####			
	#Set path where to find all .R files and folders with input/output
	statisticsfilepath <- paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/")
	model_name <- "22_hFACSELISAS_Lung_C_per_LOG_x1c_conc_filtDD_0quarterofmin_replaced"
	
	input_file <- paste0(statisticsfilepath,"input data/",model_name,".xlsx")
	parameter_labels_file <- paste0(statisticsfilepath,"input data/Parameter_labels.xlsx")
	output_folder <- paste0(statisticsfilepath,"output_Fig5C/")
	
	#read in data	
  subdataset <- data.frame(read_excel(input_file,  sheet = "Sheet 1"))
	unique_sample_ID_colname <- "Unique_Sample_ID" #column with unique identifier per row in that dataset
	ann_col <- "COPD_subclass"   #colname for colors
	ann_col2 <- "pO2_mmHg" #colname for color annotation by pO2 values
	
	sample_info <- subdataset[,1:4]
	sample_names <- subdataset[,unique_sample_ID_colname] #get samples names to later reattach to results, first column is assumed to contain sample_names or other unique identifier
	rownames(sample_info) <- sample_names
	subdataset <- subdataset[,-c(1:4)]
	rownames(subdataset) <- sample_names
	
	#read in labels
  parameter_labels <- data.frame(read_excel(parameter_labels_file,  sheet = "Parameter_labels"))
  ParamID_col <- "Parameter_Name"   #name of unique ID column with parameter names in the tabsheet Parameter_Label
  rownames(parameter_labels) <- parameter_labels[,ParamID_col]


#---Prepare for heatmaps: transform, calculate dendograms: ######			    	
  #1. scale  : column z-score scaling+centering i.e. parameters
    subdataset_map_scaled <- scale(subdataset,scale = TRUE, center = TRUE)
  
  #2. Dendograms: cluster 
    #cluster parameters, which will end up as rows in the heatmap
    mat_cluster_rows <- hclust(dist(t(subdataset_map_scaled)))     #cluster parameters, which will end up as rows in the heatmap
    mat_cluster_rows <- sort_hclust(mat_cluster_rows)  #to create HM's with sorted dendograms lowest to highest and nearest to each other
    
    #cluster samples, which will end up as columns in the heatmap
    #mat_cluster_cols <- hclust(hclust_dist)
    #mat_cluster_cols <- sort_hclust(mat_cluster_cols) #for HM's with sorted dendograms
    #no clustering by sample, keep order by pO2_mmHg
    
  #3. Calculating color break points so that color scale is centered around 0 (often needed because data often asymmetrical)
    min_boundary <- min(subdataset_map_scaled,na.rm=TRUE)
    max_boundary <- max(subdataset_map_scaled,na.rm=TRUE)
    color_breaks <- c(seq(min_boundary, min_boundary/49, -min_boundary/49), 0, seq(max_boundary/49, max_boundary, max_boundary/49)) #asymmetric scale maxing out contrast for up and down
    #breaks vector must be 1 longer than palette, because colors are valid for the interval, i.e. first color in palette is for those numbers between the first two given in the color_breaks vector, to just look at the palette
    #plots color scale alone, no numbers 
  
#---Prepare for heatmaps: set sizes and colors : ######			    	
    
  sample_annotation <- setNames(data.frame(as.numeric(rownames(subdataset_map_scaled))), unique_sample_ID_colname)
  sample_annotation <- left_join(x = sample_annotation, y = sample_info[,c(unique_sample_ID_colname, ann_col, ann_col2)],
                                   by = unique_sample_ID_colname)
  rownames(sample_annotation) <- sample_annotation[,1]
  sample_annotation <- sample_annotation[,-1, drop = FALSE]
 
  #create green color gradient, map to pO2 values
  
  pO2 <- sample_annotation[!is.na(sample_annotation[, ann_col2]), ann_col2]
  pO2_color_palette <- colorRampPalette(brewer.pal(9, "Greens"))(100)
    pO2_min <- min(pO2)
    pO2_max <- max(pO2)
    scaled_pO2 <- round(((pO2 - pO2_min) / (pO2_max - pO2_min)) * 100,0)
    scaled_pO2[1]<-1 #set zero to 1 to get lightest color
    scaled_pO2[duplicated(scaled_pO2)] <- scaled_pO2[duplicated(scaled_pO2)]+1 # to break duplicated colors

  pO2_colors <- pO2_color_palette[scaled_pO2] #to add white for the last NA pO2
  names(pO2_colors) <- pO2
  #length(scaled_pO2); length(pO2_colors)
    
  annotation_colors <- list(ann_col = c('A'='#661510', 'B'='#e34e45'),  ann_col2 = pO2_colors)
  names(annotation_colors) <- c(ann_col, ann_col2)
  
  #create names and colors for parameter annotation and nice expression labels
  parameter_labels_filt <- parameter_labels %>% filter(eval(parse(text = ParamID_col)) %in% colnames(subdataset_map_scaled)) 
          #filter to only colors of parameters in this dataset
          #reminder: must be based on subdataset_map_scaled because of possible reorders

  #grab parameter expressions, these generate nice parameter labels with super/subscrip and greek letters
  map_parameter_labels <- parameter_labels_filt[colnames(subdataset_map_scaled),"expression"]
  
  color_palette <- c(rev(colorRampPalette(brewer.pal(9, "Blues"))(48)), "#ffffff", "#ffffff" , 
                     colorRampPalette(brewer.pal(9, "Reds"))(48))	
  na_col_color = "light grey"
  border_color_set <- "dark grey"   #color of border color for every  single heatmap field, set to FALSE if none wished
  
  cell_w_spec = 8; cell_h_spec = 8; angle_col_degree= "90"
  fontsize = 8; fontsize_col = 4; fontsize_row = 6
  n_cutree_rows = 5; n_cutree_cols =  1
  
  map_w <- (nrow(subdataset_map_scaled)*cell_w_spec+350)/2.835      #is size in mm of final map
  map_h <- (nrow(t(subdataset_map_scaled))*cell_h_spec+150)/2.835	
  
  heatmap <- pheatmap(mat = t(subdataset_map_scaled), scale = "none", 
                      cluster_row = mat_cluster_rows, cluster_cols = FALSE,		#for clustering of samples give cluster object,  switch off to keep samples as sorted without clustering
                      treeheight_row = 10, treeheight_col = 15, cutree_rows = n_cutree_rows, cutree_cols =  n_cutree_cols,
                      cellheight = cell_h_spec, cellwidth = cell_w_spec,
                      annotation_col = sample_annotation,  annotation_colors = annotation_colors, annotation_row = NA,
                      color = color_palette, breaks = color_breaks,
                      fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, 
                      show_colnames = TRUE, angle_col = angle_col_degree,
                      labels_row =  parse(text=map_parameter_labels),    
                      border_color = border_color_set, na_col = na_col_color, drop_levels = TRUE) 
  ggsave(paste0(output_folder, model_name,"_heatmap.png"),plot = heatmap, device = png(), 
         width= map_w, height= map_h, units = "mm", dpi = 600);   dev.off()
  ggsave(paste0(output_folder, model_name,"_heatmap.pdf"),plot = heatmap, device = pdf(), 
         width= map_w, height= map_h, units = "mm", dpi = 600);   dev.off()
  
  