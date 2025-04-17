#calculates a random forest and exports plots of results
#---Load libraries:########--------#####
library(stringr)         #str_locate, str_sub for ParamID grabbing
library(readxl)          #read in excel file 
library(RColorBrewer)    #to add legend in MDS plots
library(ggplot2)         #for plotting 
library(ggpubr)          #ggarange, also for export
library(colorspace)      #darken/lighten
library(pheatmap)        #RF fitting heatmaps of RF
library(tidyr)           #pivot_longer
library(dplyr)           #for piping needed in functions from RF
library(caret)           #for balanced trian splits
library(randomForestSRC) #for imputation of missing values
library(randomForest)    #RF itself
library(randomForestExplainer)  #variable importance plots https://cran.rstudio.com/web/packages/randomForestExplainer/vignettes/randomForestExplainer.html

#adaptation of code from randomForestExplainer to adapt plot of min depth
# Count the trees in which each variable had a given minimal depth
min_depth_count <- function(min_depth_frame){
  tree <- NULL; minimal_depth <- NULL; variable <- NULL
  mean_tree_depth <- dplyr::group_by(min_depth_frame, tree) %>% dplyr::summarize(depth = max(minimal_depth) + 1) %>% as.data.frame()
  mean_tree_depth <- mean(mean_tree_depth$depth)
  min_depth_count <- dplyr::group_by(min_depth_frame, variable, minimal_depth) %>% dplyr::summarize(count = n()) %>% as.data.frame()
  occurrences <- stats::aggregate(count ~ variable, data = min_depth_count, sum)
  colnames(occurrences)[2] <- "no_of_occurrences"
  min_depth_count <- data.frame(variable = occurrences$variable, minimal_depth = NA,
                                count = max(min_depth_frame$tree) - occurrences$no_of_occurrences) %>% rbind(min_depth_count)
  min_depth_count <- min_depth_count[order(min_depth_count$variable, min_depth_count$minimal_depth),]
  rownames(min_depth_count) <- 1:nrow(min_depth_count)
  return(list(min_depth_count, occurrences, mean_tree_depth))
}
# Get a data frame with means of minimal depth calculated using sample = c("all_trees", "top_trees", "relevant_trees")
get_min_depth_means <- function(min_depth_frame, min_depth_count_list, mean_sample){
  .SD <- NULL; variable <- NULL
  if(mean_sample == "all_trees"){
    min_depth_count_list[[1]][is.na(min_depth_count_list[[1]]$minimal_depth), "minimal_depth"] <- min_depth_count_list[[3]]
    min_depth_means <-
      data.table::as.data.table(min_depth_count_list[[1]])[, stats::weighted.mean(.SD[["minimal_depth"]], .SD[["count"]]),
                                                           by = variable] %>% as.data.frame()
  } else if(mean_sample == "top_trees"){
    min_depth_count_list[[1]][is.na(min_depth_count_list[[1]]$minimal_depth), "count"] <-
      min_depth_count_list[[1]][is.na(min_depth_count_list[[1]]$minimal_depth), "count"] -
      min(min_depth_count_list[[1]][is.na(min_depth_count_list[[1]]$minimal_depth), "count"])
    min_depth_count_list[[1]][is.na(min_depth_count_list[[1]]$minimal_depth), "minimal_depth"] <- min_depth_count_list[[3]]
    min_depth_means <-
      data.table::as.data.table(min_depth_count_list[[1]])[, stats::weighted.mean(.SD[["minimal_depth"]], .SD[["count"]]),
                                                           by = variable] %>% as.data.frame()
  } else if(mean_sample == "relevant_trees"){
    min_depth_means <- stats::aggregate(minimal_depth ~ variable, data = min_depth_frame, mean)
  }
  colnames(min_depth_means)[2] <- "mean_minimal_depth"
  return(min_depth_means)
}

p <- function(x) {parse(text = x)} #for shorter later code cause ggplot only takes eval(parse(text= "stuff read from Excel"))

theme_set(theme_classic(base_size=8))
theme_journal <- theme(axis.line.x = element_line(color="black", linewidth = 0.5),
                       axis.line.y = element_line(color="black", linewidth = 0.5)+
                         theme(text = element_text(size=10, colour = "black")))+
  theme(axis.text.x = element_text(colour = "black"))+
  theme(axis.text.y = element_text(colour = "black"))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size = rel(1)))+ 
  theme(axis.title.y=element_text(margin=ggplot2::margin(0,5,0,0))) +  #must be called from ggplot2 because margin function masked by randomForestExplainer package
  theme(panel.grid.major = element_line(color = "light grey", linewidth=0.25, linetype = "dashed"))


#---Initializing: read in and prepare data: #####			

  set.seed(4379)   #otherwise hard to compare results due to method inherent randomness 

  #Set path where to find all .R files and folders with input/output
	statisticsfilepath <- paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/")
	model_name <- "03_hFACS_Donor_COPD_Lung_C_per_LOG_x1c_filtDD"
	input_file <- paste0(statisticsfilepath,"input data/",model_name,".xlsx")
	parameter_labels_file <- paste0(statisticsfilepath,"input data/Parameter_labels.xlsx")
	output_folder <- paste0(statisticsfilepath,"output_Fig1E_F_G/")
	
	#read in data	
  subdataset <- data.frame(read_excel(input_file,  sheet = "Sheet 1"))
	unique_sample_ID_colname <- "Unique_Sample_ID" #column with unique identifier per row in that dataset
	sample_info <- subdataset[,1:3]
	sample_names <- subdataset[,unique_sample_ID_colname] #get samples names to later reattach to results, first column is assumed to contain sample_names or other unique identifier
	rownames(sample_info) <- sample_names
	rownames(subdataset) <- sample_names
	
	#read in labels
	ParamID <- colnames(subdataset)[-c(1:3)]
  parameter_labels <- data.frame(read_excel(parameter_labels_file,  sheet = "Parameter_labels"))
  ParamID_col <- "Parameter_Name"   #name of unique ID column with parameter names in the tabsheet Parameter_Label
  rownames(parameter_labels) <- parameter_labels[,ParamID_col]
  parameter_labels <- parameter_labels[ParamID,]
  
  #get RF classification column, show distribution of samples
  ann_col <- "Diagnosis"   #colname for colors
  subdataset[,ann_col] <- as.factor(subdataset[,ann_col])         #RF need factor (cryptic error), needed for by group stratified data split
  subdataset_group <- subdataset[,ann_col]
  subdataset_noimpute <- subdataset
  #needed for RF so that it uses all Parameters
  formula_all <- as.formula(paste0(ann_col," ~ ", paste0(ParamID,collapse=" + ")))

#---Fit RF after imputing missings and splitting train data: ######			    	
  
  #if there are missings than impute them with mRF algorithm found by Tang 2017 to be best (especially well if <50% missings, missing at random and correlation existed within data)
  if(sum(is.na(subdataset[,ParamID]))!=0  ) {
    subdataset[,ParamID] <- impute.rfsrc(data = subdataset[,ParamID], mf.q = 1, verbose = TRUE) # mf.q = 0.01
  }

  train_percentage <- 0.65    #% of data split into trainingsset (for small dataset keep number around 50% otherwise too few samples in test left, for larger also 80% possible (to increase model robustness))
  intrain <- createDataPartition(y = subdataset[,ann_col], p = train_percentage, list = FALSE)
  #test <- cbind(test, intrain)
  subdataset_train <- subdataset[intrain,]
  subdataset_test <- subdataset[-intrain,]
  
  RFntree <- 5000      #number of trees to use for RF, recommend >1000 (less are often quite unstable)
  
  #DETERMINE: number of needed paramters (is done on all data), note results also change with each run, play around with stepfactor and improve cut-off for completely new datasets
  tuneRF_res <- tuneRF(x= subdataset[,ParamID], y = subdataset[,ann_col], 
                       stepFactor=1.5, plot = FALSE, ntreeTry = RFntree, trace = FALSE, improve = 0.0001)
  #stepfactor: at each iteration, mtry is inflated (or deflated) by this value; improve how much relative change is needed to continue search
  #option doBest = TRUE runs automatically random forest with optimal myTry, but with standard 500 trees
  opt_myTry <- tuneRF_res[which.min(tuneRF_res[,2]),1]  #opt_myTry <- 4
  
  ###  Plotting for RF model based on all data, for accuracy values train/test with train_percentage (e.g. 65%) train data
  RF_model <- randomForest(formula_all, data=subdataset, ntree=RFntree, mtry=opt_myTry, 
                           importance=TRUE, localImp = TRUE, proximity=TRUE) #on all data 
  RF_model_train <- randomForest(formula_all, data=subdataset_train, ntree=RFntree, 
                                 mtry=opt_myTry, importance=TRUE, localImp = TRUE, proximity=TRUE) 
        #on train data to allow independent performance evaluation with test data
  
  RF_predict <- predict(RF_model_train, subdataset_test)  #here type class also works but deliver same output, vector with predicted class
  RF_confmatr <- confusionMatrix(RF_predict, subdataset_test[,ann_col])
  
#---RF model MDS export results   ########--------#####	     
  
  plot_title <- expression("Percentage CD45"^"+"*" cells")
  
  w = 49; h = 52;  #size for nice publication figures in FACS COPD paper
  
  colors_fill <- c('COPD'='#791812', 'Donor'='#757b87')
  colors_lines <- darken(colors_fill, 0) #to create darker outer lines of shapes set <0 but <1
  names(colors_lines) <- names(colors_fill)
  
  #MDS directly from RF_model object
  RFMDS <- MDSplot(RF_model, subdataset[,ann_col], pch = as.numeric(unique(subdataset_group))) 
    #delivers non-sensical warning: In RColorBrewer::brewer.pal(nlevs, "Set1") :minimal value for n is 3, 
    #returning requested palette with 3 different levels 
  RFMDS_values <- RFMDS$points
  mds_var_per<- round(RFMDS$eig/sum(RFMDS$eig)*100,1) 
  
  plotdata_RFMDS <- data.frame(as.numeric(rownames(RFMDS_values)), MDS1 = RFMDS_values[,1], MDS2 = RFMDS_values[,2])
  colnames(plotdata_RFMDS)[1] <- unique_sample_ID_colname
  plotdata_RFMDS <- left_join(x= sample_info, y = plotdata_RFMDS, by= unique_sample_ID_colname)
#  lvls <- as.vector(unique(plotdata_RFMDS[,ann_col]))
  
  RF_Acc <- paste0("Accuracy: ",round(RF_confmatr$overall[1],2)," [",round(RF_confmatr$overall[3],2),"-",round(RF_confmatr$overall[4],2),"]")
  RFPlot_MDS <- ggplot(plotdata_RFMDS, aes(x=MDS1,y=MDS2, color = eval(p(ann_col)), fill = eval(p(ann_col)), 
                                           size = eval(p(ann_col)), shape = eval(p(ann_col)) ) ) +
    geom_point() + scale_color_manual(values=colors_lines) + scale_fill_manual(values=colors_fill) +
    scale_shape_manual(values=c(19,19)) + scale_size_manual(values=c(1,1)) + 
    xlab(paste0("MDS1 ",mds_var_per[1],"%")) + ylab(paste0("MDS2 ",mds_var_per[2],"%")) +
    annotate("text", x=min(plotdata_RFMDS$MDS1), y=max(plotdata_RFMDS$MDS2)*1.2, label = RF_Acc, size=2, hjust  ="left") +
    ggtitle(p(plot_title)) + theme_journal +
    labs(caption = paste0(model_name,"__MDS_Scores.png")) 
  
  ggsave(paste0(output_folder, "Fig1F__", model_name,"__MDS_Scores.png") ,plot = RFPlot_MDS, 
         device = png(), width=w, height=h, units = "mm", dpi = 600); dev.off()
  ggsave(paste0(output_folder, "Fig1F__", model_name,"__MDS_Scores.pdf") ,plot = RFPlot_MDS, 
         device = pdf(), width=w, height=h, units = "mm", dpi = 600); dev.off()
  
#---RF model VarImpDepth calculation and results export     ########--------#####	     	       	      	  
  
  w_vimd = 75 ; #height is set later by number of parameters
  h_vimd <- 14 + length(ParamID)*3.3
  
  #Nice minimal depth plot of all parameters
  ###caution takes 1 min to calculate the min depth distribution
  min_depth_frame <- min_depth_distribution(RF_model)   #for later plot_min_depth_distribution()
  RFplot_min_depth_all <- plot_min_depth_distribution(min_depth_frame, mean_sample = "relevant_trees", k = length(ParamID))
  
  #adaptation of code from randomForestExplainer to adapt plot of min depth
  mean_sample = "relevant_trees"; min_no_of_trees = 0; mean_round = 1
  k = length(ParamID)
  min_depth_count_list <- min_depth_count(min_depth_frame)
  min_depth_means <- get_min_depth_means(min_depth_frame, min_depth_count_list, mean_sample)
  frame_with_means <- merge(min_depth_count_list[[1]], min_depth_means)
  frame_with_means[is.na(frame_with_means$minimal_depth), "count"] <-
    frame_with_means[is.na(frame_with_means$minimal_depth), "count"] -
    min(frame_with_means[is.na(frame_with_means$minimal_depth), "count"])
  frame_with_means$mean_minimal_depth_label <-
    (frame_with_means$mean_minimal_depth - min(frame_with_means$mean_minimal_depth))/
    (max(frame_with_means$mean_minimal_depth) - min(frame_with_means$mean_minimal_depth)) *
    max(min_depth_count_list[[2]]$no_of_occurrences)
  variables <- min_depth_count_list[[2]][min_depth_count_list[[2]]$no_of_occurrences >= min_no_of_trees, "variable"]
  frame_with_means <- frame_with_means[frame_with_means$variable %in% variables, ]
  frame_with_means <- 
    within(frame_with_means, variable <-
             factor(variable, levels = unique(frame_with_means[order(frame_with_means$mean_minimal_depth), "variable"])))
  dataRF_mindepth <- frame_with_means[frame_with_means$variable %in% levels(frame_with_means$variable)[1:min(k, length(unique(frame_with_means$variable)))], ]
  dataRF_mindepth$variable <- droplevels(dataRF_mindepth$variable)
  data_for_labels <- unique(dataRF_mindepth[, c("variable", "mean_minimal_depth", "mean_minimal_depth_label")])
  data_for_labels$mean_minimal_depth <- round(data_for_labels$mean_minimal_depth, digits = mean_round)
  
  #to make Parameter_labels with super/subscripts, first pick only names left in  map
  
  RF_parameter_labels <- left_join(x = setNames(data_for_labels[,1, drop = FALSE], ParamID_col), y = parameter_labels[,c(ParamID_col,"expression")], by = ParamID_col)
  Paramter_order <- rev(levels(dataRF_mindepth$variable))
  
  RFplot_min_depth_all_png <- ggplot(dataRF_mindepth, aes(x = variable, y = count)) +
    geom_col(position = position_stack(reverse = TRUE), aes(fill = as.factor(minimal_depth))) + coord_flip() +
    scale_x_discrete(limits = Paramter_order, labels=p(RF_parameter_labels[match(Paramter_order,RF_parameter_labels[,1]),2]) ) +
    geom_errorbar(aes(ymin = mean_minimal_depth_label, ymax = mean_minimal_depth_label), size = 1.5) +
    ylab("Number of trees") + guides(fill = guide_legend(title = "")) + ggtitle(p(plot_title)) +
    geom_label(data = data_for_labels, aes(y = mean_minimal_depth_label, label = mean_minimal_depth), size = 1.5, 
               label.padding = unit(0.07, "lines")) +
    #guides(colour = guide_legend(override.aes = list(size=1))) +
    theme_journal +theme(panel.grid.major = element_line(color = "light grey", linewidth=0.25, linetype = "dashed"), 
                         axis.title.y=element_blank(),
                         legend.position="right", legend.key.size = unit(0.25,"line"), 
                         axis.text.x = element_text(size = 5),  axis.text.y = element_text(size = 6, color ="black") ) +
    labs(caption = paste0(model_name,"__VarImpminDepth.png")) 
  
  # legend.position="bottom", legend.box = "vertical", legend.key.size = unit(0.25,"line"))
  
  ggsave(paste0(output_folder, "Fig1E__", model_name,"__VarImpminDepth.png"), 
         plot = RFplot_min_depth_all_png, device = png(), width=w_vimd, height=h_vimd, units = "mm", dpi = 600)
  dev.off()
  ggsave(paste0(output_folder, "Fig1E__", model_name,"__VarImpminDepth.pdf"), 
         plot = RFplot_min_depth_all_png, device = pdf(), width=w_vimd, height=h_vimd, units = "mm", dpi = 600)
  dev.off()

#---Plot Heatmap of log fold changes for each parameter (right side of VarImpDepth plot)     ########--------#####	     	       	      	  
  
  #get log fold changes (LFC), for current code simplicity only second lvl vs. first lvl
  group_medians <- setNames(aggregate(10^subdataset_noimpute[,ParamID], 
                                      list(group = subdataset_noimpute[,ann_col]), 
                                      FUN = function(x) median(x, na.rm = TRUE)), #have to set na.rm specifically
                            c(ann_col, ParamID) )
  
  LFC_medians <- data.frame(colnames(group_medians[,-1]), t(group_medians[,-1]))
  colnames(LFC_medians) <- c(ParamID_col, as.character(group_medians[,1]))
  ref_group <- "Donor"
  lvls_to_calc <- lvls[-which(lvls == ref_group)]  #calculate LFC only for non-reference group
  n_lvls <- length(lvls)
  
  #create additional columns with same parameter order as in Min_depth and merge as y (=position on y axis)
  temp_min_depth_for_sort <- unique(dataRF_mindepth[,c("variable","mean_minimal_depth")])
  temp_min_depth_for_sort <- temp_min_depth_for_sort[order(temp_min_depth_for_sort$mean_minimal_depth),]
  y_merge <- data.frame(Parameter_Name = as.character(temp_min_depth_for_sort[,1]), y=nrow(temp_min_depth_for_sort):1)
  y_labels <- as.character(y_merge[,1])
  names(y_labels) <- as.character(y_merge[,2])
  y_labels <- RF_parameter_labels[match(y_labels,RF_parameter_labels[,1]),2]
  
  LFC_medians_long  <- data.frame(LFC_medians, comparison = paste0(lvls_to_calc[1]," vs ",ref_group), 
                                  FC = LFC_medians[,lvls_to_calc[1]]/LFC_medians[,ref_group],
                                  LFC = log(LFC_medians[,lvls_to_calc[1]]/LFC_medians[,ref_group],2)) #log2 of ratio group/ref_group
  LFC_medians_long <- data.frame(LFC_medians_long, comp_group_higher = cut(LFC_medians_long[,"LFC"], 
                                                                           c(-Inf,0,Inf), 
                                                                           c(ref_group,lvls_to_calc[1])) )#assign in additional column which group had higher levels (for coloring by group)
  
  LFC_medians_long <- left_join(x = LFC_medians_long, y = y_merge, by = ParamID_col)  #match back y positions of parameters (to have same order as in RF var Imp plot)
  
  RFplot_groupLFC_map <- 
    ggplot(LFC_medians_long, aes(x = comparison, y = y)) +
    geom_tile(aes(color = comp_group_higher, fill = comp_group_higher), size = 0.5, width = 0.75, height = 0.75) +
    geom_text(aes(label=round(LFC,1)), colour="white", size = 2)+
    scale_color_manual(values=colors_fill) + scale_fill_manual(values=colors_fill) + 
    scale_x_discrete(position = "top") + 
    scale_y_continuous(breaks = y_merge[,2], labels=p(y_labels)) +
    theme_void() + 
    theme(axis.text.x = element_text(colour = "black",size=6)) +  #, plot.margin = margin(0, 0, 0, 0, "cm"), angle=45,  hjust = 0, vjust = 0
    theme(axis.text.y = element_text(colour = "black",size=6, hjust = 1)) +
    guides(fill = guide_legend(title = "")) +  ggtitle(p(plot_title)) +
    theme(legend.position="none")
  
  ggsave(paste0(output_folder, "Fig1E_right_", model_name,"__groupLFC_map.png"), 
         plot = RFplot_groupLFC_map, device = png(), width=w_vimd-34+4*(n_lvls-1), height=h_vimd+4, units = "mm", dpi = 600)
  dev.off() 
  ggsave(paste0(output_folder, "Fig1E_right_",model_name,"__groupLFC_map.pdf"), 
         plot = RFplot_groupLFC_map, device = pdf(), width=w_vimd-34+4*(n_lvls-1), height=h_vimd+4, units = "mm", dpi = 600)
  dev.off() 
  
#---Fit RF on only cells in >300 tree at root node: ######			    	
  ParamID_top300T <- c("CD3", "CD4", "CD8", "CD19", "Neutrophils", "Mono_int", "Mono_clas")
  subdataset_top300T <- data.frame(sample_info, subdataset_noimpute[,ParamID_top300T])
  
  #if there are missings than impute them with mRF algorithm found by Tang 2017 to be best (especially well if <50% missings, missing at random and correlation existed within data)
  if(sum(is.na(subdataset_top300T[,ParamID_top300T]))!=0  ) {
    subdataset_top300T[,ParamID_top300T] <- impute.rfsrc(data = subdataset_top300T[,ParamID_top300T], mf.q = 1, verbose = TRUE) 
  }
  subdataset_top300T[,ann_col] <- as.factor(subdataset_top300T[,ann_col])
  formula_all_top300T <- as.formula(paste0(ann_col," ~ ", paste0(ParamID_top300T,collapse=" + ")))
  
  
  intrain_top300T <- createDataPartition(y = subdataset_top300T[,ann_col], p = train_percentage, list = FALSE)
  subdataset_train_top300T <- subdataset_top300T[intrain,]
  subdataset_test_top300T <- subdataset_top300T[-intrain,]
  
  #DETERMINE: number of needed paramters (is done on all data), note results also change with each run, play around with stepfactor and improve cut-off for completely new datasets
  tuneRF_res_top300T <- tuneRF(x= subdataset_top300T[,ParamID_top300T], y = subdataset_top300T[,ann_col], 
                       stepFactor=1.5, plot = FALSE, ntreeTry = RFntree, trace = FALSE, improve = 0.0001)
  #stepfactor: at each iteration, mtry is inflated (or deflated) by this value; improve how much relative change is needed to continue search
  #option doBest = TRUE runs automatically random forest with optimal myTry, but with standard 500 trees
  opt_myTry_top300T <- tuneRF_res_top300T[which.min(tuneRF_res_top300T[,2]),1]  #opt_myTry <- 4
  
  ###  Plotting for RF model based on all data, for accuracy values train/test with train_percentage (e.g. 65%) train data
  RF_model_top300T <- randomForest(formula_all_top300T, data=subdataset_top300T, ntree=RFntree, mtry=opt_myTry_top300T, 
                           importance=TRUE, localImp = TRUE, proximity=TRUE) #on all data 
  RF_model_train_top300T <- randomForest(formula_all_top300T, data=subdataset_train_top300T, ntree=RFntree, 
                                 mtry=opt_myTry_top300T, importance=TRUE, localImp = TRUE, proximity=TRUE) 
  #on train data to allow independent performance evaluation with test data
  
  RF_predict_top300T <- predict(RF_model_train_top300T, subdataset_test_top300T)  #here type class also works but deliver same output, vector with predicted class
  RF_confmatr_top300T <- confusionMatrix(RF_predict_top300T, subdataset_test_top300T[,ann_col])
  
#---RF model MDS export results   ########--------#####	     
  
  #MDS directly from RF_model object
  RFMDS_top300T <- MDSplot(RF_model_top300T, subdataset_top300T[,ann_col], pch = as.numeric(unique(subdataset_group))) 
  #delivers non-sensical warning: In RColorBrewer::brewer.pal(nlevs, "Set1") :minimal value for n is 3, 
  #returning requested palette with 3 different levels 
  RFMDS_values_top300T <- RFMDS_top300T$points
  mds_var_per_top300T <- round(RFMDS_top300T$eig/sum(RFMDS_top300T$eig)*100,1) 
  
  plotdata_RFMDS_top300T <- data.frame(as.numeric(rownames(RFMDS_values_top300T)), 
                                       MDS1 = RFMDS_values_top300T[,1], MDS2 = RFMDS_values_top300T[,2])
  colnames(plotdata_RFMDS_top300T)[1] <- unique_sample_ID_colname
  plotdata_RFMDS_top300T <- left_join(x= sample_info, y = plotdata_RFMDS_top300T, by= unique_sample_ID_colname)

  RF_Acc_top300T <- paste0("Accuracy: ",round(RF_confmatr_top300T$overall[1],2),
                   " [",round(RF_confmatr_top300T$overall[3],2),"-",
                   round(RF_confmatr_top300T$overall[4],2),"]")
  RFPlot_MDS_top300T <- ggplot(plotdata_RFMDS_top300T, aes(x=MDS1,y=MDS2, color = eval(p(ann_col)), fill = eval(p(ann_col)), 
                                           size = eval(p(ann_col)), shape = eval(p(ann_col)) ) ) +
    geom_point() + scale_color_manual(values=colors_lines) + scale_fill_manual(values=colors_fill) +
    scale_shape_manual(values=c(19,19)) + scale_size_manual(values=c(1,1)) + 
    xlab(paste0("MDS1 ",mds_var_per[1],"%")) + ylab(paste0("MDS2 ",mds_var_per[2],"%")) +
    annotate("text", x=min(plotdata_RFMDS_top300T$MDS1), y=max(plotdata_RFMDS_top300T$MDS2)*1.2, 
             label = RF_Acc_top300T, size=2, hjust  ="left") +
    ggtitle(p(plot_title)) + theme_journal +
    labs(caption = paste0(model_name,"__MDS_Scores_top300T.png")) 
  
  ggsave(paste0(output_folder, "Fig1G__", model_name,"__MDS_Scores_top300T.png") ,plot = RFPlot_MDS_top300T, 
         device = png(), width=w, height=h, units = "mm", dpi = 600); dev.off()
  ggsave(paste0(output_folder, "Fig1G__", model_name,"__MDS_Scores_top300T.pdf") ,plot = RFPlot_MDS_top300T, 
         device = pdf(), width=w, height=h, units = "mm", dpi = 600); dev.off()
  
  