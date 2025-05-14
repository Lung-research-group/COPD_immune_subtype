# *Machine Learning-Assisted Multimodal Analysis Identifies an Immune-Driven COPD Subtype Associated with Increased Emphysema Across GOLD Stages*

Natalie Bordag<sup>§</sup>, Katharina Jandl<sup>§</sup>, Ayu Hutami Syarif<sup>§</sup>, Jürgen Gindlhuber, Diana Schnoegl, Ayse Ceren Mutgan, Vasile Foris, Konrad Hoetzenecker, Panja M. Boehm, Robab Breyer-Kohansal, Katarina Zeder, Gregor Gorkiewicz, Francesca Polverino, Slaven Crnkovic, Grazyna Kwapiszewska, Leigh M. Marsh<sup>$</sup>

§ These authors contributed equally to the work and are listed alphabetically  
$ Corresponding author


<p align="center"><img src="/COPD_graphical_abstract.jpg" height="700" width=900"></p>

## Description
This repository contains the R codes used to generate all the figures in the article *Machine Learning-Assisted Multimodal Analysis Identifies an Immune-Driven COPD Subtype Associated with Increased Emphysema Across GOLD Stages*

The raw data used for the Flow Cytometry and Cytokine analyses have been deposited to Mendeley Data (human FACS data, accession number [XXXXX]). The raw data for scRNA-seq  were generated from control and COPD lungs, published by [Sauler et al. (2022)] (https://doi.org/10.1038/s41467-022-28062-9) (GEO accession: GSE136831). The raw data for spatial transcriptomics were produced from  GOLD1-2 and GOLD 3-4 COPD lungs, published by Polverino, et al (GEO accession: XXX). 

## Content

* `Figure 1-S2/`: Scripts and pre-processed data to reproduce Figure 1
  * `Fig1_C_D_E_F_G/`: Scripts to reproduce figures 1C-G
    * `input data/`: Input tables to load into R codes (R_Fig1C_PCA_biplot.r, R_Fig1D_OPLSDA_Scores_plot.r, R_Fig1E_F_G_RandomForest_Vardepth_MDS_plots.r)
  * `iScience_input_data/`: Input tables to load into R codes (Fig1B-1H-S2_COPD.R) 
* `Figure 2/`: Scripts and pre-processed data to reproduce Figure 2
    * `input data/`: Input tables to load into R codes (Fig2I-2J_COPD.R)
    * `iScience_input_data/`: Input tables to load into R codes (Fig2B-2C-2D-2E-2F-2G-S6-S7-S8_COPD.R) 
* `Figure 3/`: Scripts and pre-processed data to reproduce Figure 3
    * `input data/`: Input tables to load into R codes (Fig3_COPD.R)
* `Figure 4/`: Scripts and pre-processed data to reproduce Figure 4
    * `Fig4Bdown/`
      *`input data/`: Input tables to load into R codes (R_Fig4Bdown_OPLSDA_Scores_plot.r, R_Fig4_COPD_subclass.r)
    * `Fig4Cdown-4D/`
      *`input data/`: Input tables to load into R codes (R_Fig4Cdown_OPLSDA_Scores_plot.r, R_Fig4D_Heatmap.r)
    * `iScience_input_data/`: Input tables to load into R codes (Fig4B-4C top.R, Fig4E-4F-4G-5B_COPD.R) 
* `Figure 5/`: Scripts and pre-processed data to reproduce Figure 5
    * `Fig5C/`
      *`input data/`: Input tables to load into R codes (R_Fig5C_Heatmap.r)
    * `Fig5D_E/`
      *`input data/`: Input tables to load into R codes (Fig5D_network.R)
    * `iScience_input_data/`: Input tables to load into R codes (Fig4E-4F-4G-5B_COPD.R)
* `Figure 6/`: Scripts and pre-processed data to reproduce Figure 6
    *`input data/`: Input tables to load into R codes (Fig6_ST.R)
* `Figure S10/`: Scripts and pre-processed data to reproduce Figure S10
    *`iScience_input_data/`: Input tables to load into R codes (Fig_S10.R)
* `Figure S11_S12/`: Scripts and pre-processed data to reproduce Figure S11-12
    *`input data/`: Input tables to load into R codes (Fig_S11_S12.R)
* `Figure S3/`: Scripts and pre-processed data to reproduce Figure S3
    *`iScience_input_data/`: Input tables to load into R codes (FigS3_COPD.R)
* `Figure S6-S7-S8/`: Scripts and pre-processed data to reproduce Figure S3, S6, S8
    *`iScience_input_data/`: Input tables to load into R codes (Fig2-S6-S7-S8_COPD.R)
* `Input Data Files/`: all input data tables to recreate figures
    *`FACS/`: Input tables for flow cytometry analysis
    *`GeoMx/`: Input tables for spatial transcriptomics analysis
    *`scRNA/`: Input tables for scRNA-seq analysis
* `Metadata/`: Metadata of individual cohorts

  
## Get started

If you want to reproduce some of the figures presented here, it is recommended to do the following steps:

1. Move the `semla_analysis` folder and content to another location (this analysis path has its own R project file and package requirements)  
2. Download the spaceranger output data and/or STUtility/Seurat objects and place them in their correct locations (see `data/README.txt`)  
3. Install [`renv`](https://rstudio.github.io/renv/articles/renv.html) and install the necessary packages from the `renv.lock` file by using `renv::restore()`
4. Manually install the packages within the `bin/` folder(s) and [STUtility v. 1.1.1](https://github.com/jbergenstrahle/STUtility/releases/tag/1.1.1)  
5. Read the `scripts/README.txt` for information about all scripts and their content


## Contact

For further questions related to this repository, its content, and the paper, please contact Leigh M. Marsh (Leigh.Marsh@medunigraz.at) Telephone: +43 316 385 72911
