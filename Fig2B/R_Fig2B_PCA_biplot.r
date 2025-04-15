#calculates a classical PCA (prcomp) with missings imputated by EM
#---Load libraries:########--------#####
	library(readxl)     ##read in excel file 
  library(missMDA)    ##imputePCA to impute missings in dataset to enable later classical PCA
  library(ggplot2)    ##for plotting 
	library(ggpubr)     ##for plotting 
  library(colorspace) ##darken/lighten
	library(factoextra) ##for PCA loadings plots with arrows 

#---Initializing:  setting all names, sizes, colors, read in data etc.: #####			
	#Set path where to find all .R files and folders with input/output
	statisticsfilepath <- paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/")
	model_name <- "05_hELISAL_Donor_COPD_conc_tprotNorm_OEME_LOG_0quarterofmin_replaced"
	input_file <- paste0(statisticsfilepath,"input data/",model_name,".xlsx")
	parameter_labels_file <- paste0(statisticsfilepath,"input data/Parameter_labels.xlsx")
	output_folder <- paste0(statisticsfilepath,"output_Fig2B/")
	
	#read in data	
  subdataset <- data.frame(read_excel(input_file,  sheet = "Sheet 1"))
	unique_sample_ID_colname <- "Unique_Sample_ID" #column with unique identifier per row in that dataset
	sample_info <- subdataset[,1:3]
	sample_names <- subdataset[,unique_sample_ID_colname] #get samples names to later reattach to results, first column is assumed to contain sample_names or other unique identifier
	rownames(sample_info) <- sample_names
	subdataset <- subdataset[,-c(1:3)]
	rownames(subdataset) <- sample_names
	
	#read in labels
  parameter_labels <- data.frame(read_excel(parameter_labels_file,  sheet = "Parameter_labels"))
  ParamID_col <- "Parameter_Name"   #name of unique ID column with parameter names in the tabsheet Parameter_Label
  rownames(parameter_labels) <- parameter_labels[,ParamID_col]

#---Impute missings and perform classical PCA: ######
    if(sum(is.na(subdataset))!=0) { #impute if missings exist
      if(ncol(subdataset)<=6) {ncomp_impute <- ncol(subdataset)-1} else {ncomp_impute <- 6}
      subdataset_imputed <- imputePCA(subdataset,ncomp_impute)$completeObs   
    } else {
      subdataset_imputed <- subdataset
    }

    PCA_modelobject <- prcomp(subdataset_imputed, scale=TRUE, center=TRUE) #classical PCA calculation

#---Prepare results plotting: ######			    	
  #set names, shapes, sizes, colors
  plot_title <- c("Lung")
  ann_col <- "Diagnosis"   #colname for colors
  colors_fill <- c('Donor'='#757b87', 'COPD'='#791812')
  colors_lines <- darken(colors_fill, 0) #to create darker outer lines of shapes set <0 but <1
  names(colors_lines) <- names(colors_fill)

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
  
#---Export Bi-plots: ######			    	
  
  subdataset_imputed_renamed <- subdataset_imputed   #is a matrix with samples in rows and parameters in columns
  biplot_param_names <- colnames(subdataset_imputed)
  biplot_param_names <- data.frame(parameter_labels[biplot_param_names,"Plot_label"])
  colnames(subdataset_imputed_renamed) <- eval(parse(text = biplot_param_names))
  
  PCA_modelobject <- prcomp(subdataset_imputed_renamed, scale=TRUE, center=TRUE) #classical PCA calculation
  
  xlabel_PC2 <- paste0("PC2 (", round(100*summary(PCA_modelobject)$importance[2,2],1), "%)"); 
  ylabel_PC4 <- paste0("PC4 (", round(100*summary(PCA_modelobject)$importance[2,4],1), "%)");
        
  Biplot_cPCA_PC24 <- fviz_pca_biplot(PCA_modelobject,  axes = c(2, 4), title = plot_title, 
                                      habillage= sample_info[rownames(PCA_modelobject$x), ann_col],
                                       mean.point = FALSE, pointsize = 1, pointshape = 19, 
                                       label ="var", select.var = list(contrib = 10), #addEllipses=TRUE, ellipse.level=0.95
                                       labelsize = 1, arrowsize = 0.1, col.var = "#252525", repel = TRUE) + #alpha.var="contrib",  #can't color by contrib since uses already color_scale for PCA scores
          scale_color_manual(values=colors_lines) +
          expand_limits(x = c(min(PCA_modelobject$x[,"PC2"])*1.1, max(PCA_modelobject$x[,"PC2"])*1.1),
                        y = c(min(PCA_modelobject$x[,"PC4"])*1.1, max(PCA_modelobject$x[,"PC4"])*1.1)) +  
          guides(color = guide_legend(override.aes = list(size = 1))) + #reduces dot size in legend
          xlab(xlabel_PC2) + ylab(ylabel_PC4) + theme_journal +   # +coord_fixed() #makes quadratic
          theme(legend.position="none", plot.title = element_text(color="black", size = 6),
                plot.caption = element_text(size = 4, hjust = 0), plot.caption.position = "plot",
                axis.title.x = element_text(colour = "black", size = 6), 
                axis.title.y = element_text(colour = "black", size = 6) ) + 
          labs(caption = paste0(model_name,"__cPCA_biplot.png")) 
  ggsave(paste0(output_folder, model_name,"__cPCA_PC24_biplot.png"),plot = Biplot_cPCA_PC24 , device = png(), width=w+5, 
         height=h+7, units = "mm", dpi = 900); dev.off()
  ggsave(paste0(output_folder, model_name,"__cPCA_PC24_biplot.pdf"),plot = Biplot_cPCA_PC24 , device = pdf(), width=w+5, 
         height=h+7, units = "mm", dpi = 900); dev.off()
        