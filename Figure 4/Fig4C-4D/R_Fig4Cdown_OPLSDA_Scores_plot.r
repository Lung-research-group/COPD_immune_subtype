#calculates a OPLS-DA from the MeaboAnalystR package
#---Load libraries:########--------#####
	library(readxl)     ##read in excel file 
  library(missMDA)    ##imputePCA to impute missings in dataset to enable later classical PCA
  library(ggplot2)    ##for plotting 
	library(ggpubr)     ##for plotting 
  library(colorspace) ##darken/lighten
  library(MetaboAnalystR)  #for OPLS-DA

#---Initializing:  setting all names, sizes, colors, read in data etc.: #####			
	#Set path where to find all .R files and folders with input/output
	statisticsfilepath <- paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/")
	model_name <- "22_hFACSELISAS_Lung_C_per_LOG_x1c_conc_filtDD_0quarterofmin_replaced"
	
	input_file <- paste0(statisticsfilepath,"input data/",model_name,".xlsx")
	parameter_labels_file <- paste0(statisticsfilepath,"input data/Parameter_labels.xlsx")
	output_folder <- paste0(statisticsfilepath,"output_Fig4Cdown/")
	  if (!dir.exists(output_folder)){dir.create(output_folder)}
	output_folder_MetaboAnalystR <- paste0(output_folder,"MetaboAnalystR_files/")
	if (!dir.exists(output_folder_MetaboAnalystR)){dir.create(output_folder_MetaboAnalystR)}
	setwd(output_folder_MetaboAnalystR) #needed for MetaboAnalystR package to work
	
	#read in data	
  subdataset <- data.frame(read_excel(input_file,  sheet = "Sheet 1"))
	unique_sample_ID_colname <- "Unique_Sample_ID" #column with unique identifier per row in that dataset
	ann_col <- "COPD_subclass"   #colname for colors/group calculation
	
	sample_info <- subdataset[,1:3]
	sample_names <- subdataset[,unique_sample_ID_colname] #get samples names to later reattach to results, first column is assumed to contain sample_names or other unique identifier
	rownames(sample_info) <- sample_names
	subdataset <- subdataset[,-c(1:3)]
	rownames(subdataset) <- sample_names
	subdataset <- as.matrix(subdataset)   #matrix needed for later calculations,
	subdataset_group <- sample_info[,ann_col]  
	
	#read in labels
  parameter_labels <- data.frame(read_excel(parameter_labels_file,  sheet = "Parameter_labels"))
  ParamID_col <- "Parameter_Name"   #name of unique ID column with parameter names in the tabsheet Parameter_Label
  rownames(parameter_labels) <- parameter_labels[,ParamID_col]

#---Calculate MetaboAnalystR OPLS-DA models ########--------#####	  
  
  #initialize MetaboAnalystR object to contain data and results
  mSet <- InitDataObjects("pktable", "stat", FALSE)   
  
  #assemble data for specific model, safe as .csv, reload into MetaboanalystR package
  MetaboAnalyst_data <- data.frame(sample_names, subdataset_group, subdataset)
  colnames(MetaboAnalyst_data)[2] <- ann_col
  write.csv(MetaboAnalyst_data, paste0(output_folder_MetaboAnalystR,"MetaboAnalyst_data_",model_name,".csv"), row.names = FALSE)
  mSet <- Read.TextData(mSet, paste0(output_folder_MetaboAnalystR,"MetaboAnalyst_data_",model_name,".csv"), "rowu", "disc")   
    #reads into View(mSet$dataSet$orig)
  mSet <- SanityCheckData(mSet)                     #uses View(mSet$dataSet$orig) to create View(mSet$dataSet$preproc)
  mSet <- RemoveMissingPercent(mSet, percent=0.5)   #uses View(mSet$dataSet$preproc), writes back to View(mSet$dataSet$preproc)
  mSet <- ImputeMissingVar(mSet, method="knn_var")    #unlcear why sometimes one or rather line works
  
  #make data-transformations, here only Autoscaling for later OPLS-DA needed, other data-transformation like Logx+1 were done before
  mSet<-PreparePrenormData(mSet)
  mSet<-Normalization(mSet, "NULL", "NULL", "AutoNorm", ratio=FALSE, ratioNum=20)   #writes mSetObj$dataSet$proc (im?putes based on mSetObj$dataSet$preproc)
  
  scaling <- "AutoNorm"
  #OPLS-DA itself
  mSet<-OPLSR.Anal(mSet, reg=TRUE) #reg Add reg (regression i.e. if class order matters), creates 2 Excels in R home directory
  mSet<-PlotOPLS2DScore(mSet, paste0(model_name, "_", scaling, "_opls_score2d_"), format = "png", dpi=600, width=NA, 1,2,0.95,1,0)    #Create a 2D oPLS-DA score plot
  mSet<-PlotOPLS.Splot(mSet, paste0(model_name, "_", scaling, "_opls_splot_"), format = "png", dpi=600, width=NA)     # Create a significant features plot
  mSet<-PlotOPLS.MDL(mSet, paste0(model_name, "_", scaling, "_opls_mdl_"), format = "png", dpi=600, width=NA)   # Create a plot of the model overview
  mSet$analSet$opls.reg <- TRUE   #occurred as error message since update to R version 4
  mSet<-OPLSDA.Permut(mSet, 1000)    # Perform oPLS-DA permutation, small dataset 1000 is fine, larger better 100 takes too long otherwise
  mSet<-PlotOPLS.Permutation(mSet, paste0(model_name, "_", scaling, "_opls_perm_"), format = "png", dpi=600, width=NA)    # plot oPLS-DA permutation 
  
  
#---Prepare results plotting: ######			    	
  #set names, shapes, sizes, colors
  plot_title <- expression("Percentage CD45"^"+"*" cells, Plasma")
  colors_fill <- c('A'='#661510', 'B'='#e34e45')
  colors_lines <- darken(colors_fill, 0) #to create darker outer lines of shapes set <0 but <1
  names(colors_lines) <- names(colors_fill)

  theme_set(theme_pubr(base_size=6, border = T, legend = "none"))
  theme_journal <- theme(panel.grid.major = element_line(color = "light grey", linewidth = 0.2, linetype = "solid"), 
                         panel.grid.minor = element_line(color = "light grey", linewidth = 0.05, linetype = "solid"), 
                         plot.title = element_text(size = rel(1)),
                         axis.line.x = element_line(color="black", linewidth = 0.2),
                         axis.line.y = element_line(color="black", linewidth = 0.2),
                         axis.text = element_text(color="#393939", size = 6),
                         panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.2)
  )
 
  w = 49; h = 52;  #size settings for plot exports

#---: ######			    	
  
  xlabel <- paste0("T score 1 (", round(100*mSet$analSet$oplsda$modelDF["p1", "R2X"],1), "%)");
  ylabel <- paste0("Orth T score 1 (", round(100*mSet$analSet$oplsda$modelDF["o1", "R2X"],1), "%)");
  R2Y <- mSet$analSet$oplsda$perm.res$r.vec[1]
  Q2 <- mSet$analSet$oplsda$perm.res$q.vec[1]
  
  df_temp_merge <- data.frame(rownames(mSet$analSet$oplsda$scoreMN), mSet$analSet$oplsda$scoreMN)
  colnames(df_temp_merge)[1] <- unique_sample_ID_colname
  plotdata_OPLS <- merge(x= sample_info, y = df_temp_merge, by= unique_sample_ID_colname)
  df_temp_merge <- data.frame(rownames(mSet$analSet$oplsda$orthoScoreMN), mSet$analSet$oplsda$orthoScoreMN)
  colnames(df_temp_merge)[1] <- unique_sample_ID_colname
  plotdata_OPLS <- merge(x= plotdata_OPLS, y = df_temp_merge, by= unique_sample_ID_colname)
  
  plotdata_OPLS[,ann_col] <- as.factor(plotdata_OPLS[,ann_col])
    
  #calculate x/y axis limits by also calculating the ellipse for the current coloring groups ellipse, 
  #caution ellipses are drawn to coloring, NOT to the discriminator, Code directly copied/adapted from MetaboAnalystR function PlotOPLS2DScore
    lvs <- unique(plotdata_OPLS[,ann_col]) 
    pts.array <- array(0, dim=c(100,2,length(lvs)));
    for(i in 1:length(lvs)){
      inx <- plotdata_OPLS[,ann_col] == lvs[i];
      groupVar <- var(cbind(plotdata_OPLS[inx,"p1"],plotdata_OPLS[inx,"o1"]), na.rm=T);
      groupMean <- cbind(mean(plotdata_OPLS[inx,"p1"], na.rm=T),mean(plotdata_OPLS[inx,"o1"], na.rm=T));
      pts.array[,,i] <- ellipse::ellipse(groupVar, centre = groupMean, level = 0.95, npoints=100);
    }
    
    xrg <- range(plotdata_OPLS[,"p1"], pts.array[,1,]);
    yrg <- range(plotdata_OPLS[,"o1"], pts.array[,2,]);
    x.ext<-(xrg[2]-xrg[1])/12;
    y.ext<-(yrg[2]-yrg[1])/12;
    xlims<-c(xrg[1]-x.ext, xrg[2]+x.ext);
    ylims<-c(yrg[1]-y.ext, yrg[2]+y.ext);
    
    
    scoresPlot_OPLS <- ggplot(data = plotdata_OPLS, aes(x=p1,y=o1, 
                                                        color = eval(parse(text = ann_col)), 
                                                        fill = eval(parse(text = ann_col)), 
                                                        size = eval(parse(text = ann_col)), 
                                                        shape = eval(parse(text = ann_col)) ) ) +
      geom_point() + 
      stat_ellipse(aes(group = eval(parse(text = ann_col))), geom="polygon", alpha = 0.3, linewidth = 0.25) +
      scale_color_manual(values=colors_lines) + scale_fill_manual(values=colors_fill) +
      scale_shape_manual(values=c(19,19)) + scale_size_manual(values=c(1,1)) + 
      xlab(xlabel) + ylab(ylabel) +expand_limits(x = xlims, y= ylims) + ggtitle(parse(text = plot_title)) + 
      theme_journal +
      ggplot2::annotate("text", x=xlims[1]*0.75, y=ylims[2], label = bquote(R^2 * 'Y' == .(R2Y * 100) * '%'), size = 2)+  
      ggplot2::annotate("text", x=xlims[1]*0.78, y=ylims[2]*0.85, label = bquote(Q^2 * 'Y' == .(Q2 * 100) * '%'), size = 2)  
  
  suppressWarnings(ggsave(paste0(output_folder, model_name,"__OPLS_Scores___",ann_col,".png"),
                          plot = scoresPlot_OPLS, device = png(), width=w, height=h, units = "mm", dpi = 600))
  dev.off()
  #and safe a .pdf of that
  suppressWarnings(ggsave(paste0(output_folder, model_name,"__OPLS_Scores___",ann_col,".pdf"),
                          plot = scoresPlot_OPLS, device = pdf(), width=w, height=h, units = "mm", dpi = 600))
  dev.off()
  
  