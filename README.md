# DARs

Code and data for the paper:

Matthews, T.J. et al. (2023) A global analysis of avian island diversity–area relationships in the Anthropocene. Ecology Letters, in press.

Code authors: Tom Matthews and Francois Rigal (Figure 5 code)

## Software used

Analyses run using the University of Birmingham's BLUE BEAR HPC,
using 51 cores and 500gb memory.

R version 4.2.0 (2022-04-22)

R packages:

ape_5.6-2  
picante_1.8.2 (github repo: "txm676/picante@98fdbde"; archived DOI: https://doi.org/10.5281/zenodo.7596317)    
dplyr_1.0.8   
sars_1.3.6  
foreach_1.5.2  
doParallel_1.0.17  
cluster_2.1.3  
AICcmodavg_2.3-1  
VGAM_1.1-6  
ggplot2_3.3.5  
nlme_3.1-157  
glmmTMB_1.1.3   
car_3.0-12    
tidyr_1.2.0  
colorBlindness_0.1.9  
reshape2_1.4.4  
ppcor_1.1  
lemon_0.4.5


## Restrictions and embargos
1)	In the Wakatobi dataset (O’Connell et al. 2020), we have removed one species from the published dataset to avoid publishing the islands on which it is found: Cacatua sulphurea – a critically endangered cockatoo that is at risk from the caged bird trade. Please contact Tom Martin (a paper co-author; tom_martin_2010@yahoo.co.uk) if this information is required.

2)	The Zhoushan dataset (Zhao et al.) has only recently been developed through an extensive field work campaign. As such, this dataset has been embargoed to allow for it to be used in a different publication. The dataset has been archived at Zenodo with an embargo period of 12 months: 10.5281/zenodo.7596594 (https://zenodo.org/record/7596594)

## R code files
Within the R Code directory, there are four .R files:

1)	DAR_modelling – this is the main script for running the analyses, secondary analyses and making the figures. It uses parallel processing. Note here, the arguments in the main fit_DARs function have been set to allow the scripts to be run through v quickly but these are NOT the arguments used in the main analyses. To replicate these, set grid_start = "exhaustive", grid_n = 25000, null_n = 999, prune_trees = FALSE. But beware this takes many hours (days) to run!

2)	DAR_SOURCE – this holds all of the functions that are used in DAR_modelling.

3)	DAR_DATA_SOURCE – loads in the consensus phylogeny, the trait data (for extant and extinct species) and the datasets. Depending on what settings are used in DAR_modelling, uses body-size corrected traits or not, and the full datasets or the land bird versions.

4)	HEATMAP – the code for making Figure 5.

## Data files
A separate Dataset_information document is included as part of the repository which provides information for each dataset individually.
Within the Data directory, there are three folders:

1)	Island_datasets – the full datasets used in the main analyses. Datasets are stored as .csvs. The format of the datasets is a presence-absence matrix, with species as rows and sites as columns (this is transposed in the scripts). The bottom two rows include the island area and species richness (sp.r) totals for each island. All island areas are in hectares (which are converted to km2 in the scripts).  All species names follow the Jetz et al (2012) BirdTree taxonomy. The first folder contains the 25 habitat island datasets, with the third folder representing the true island datasets. The second folder (Land_bird_datasets) contains all the dataset versions that only include land birds (excluding the the marine, coastal, wetland and riverine species).

Within the True_island_datasets folder, there is an "Alternative_versions" folder which contains the dataset versions with alien species included, and those representing the historic and prehistoric periods (following the same format as above).

2) Predictors - a single csv file that contains the elevation (Elev_max), climate (Bio1_m, Bio12_m), isolation (Iso), MeanDist and island type (Type_) information for all datasets.

3) Species_datasets -	contains two csv files and the consensus phylogeny:

•	The trait data for extant species (BirdTree taxonomy). Stored as a csv with 9,993 rows and 12 columns (the 9 traits in the analyses and family, order and habitat data). Rows are species and columns trait information. Data are from AVONET (https://onlinelibrary.wiley.com/doi/full/10.1111/ele.13898) and users should cite this paper if using these data in any publication (it is probably preferable to source the data direct from AVONET in this case). See AVONET metadata for information on all traits.

•	The trait data for extinct species. A csv with 159 rows and 15 columns. Rows are species and columns trait information. Contains data for the same 9 traits as for the extant species, in addition to the Archipelago the species was found on, the Order, Family and Genus, as well as whether or not the species is Aquatic (TRUE = aquatic). Note that some trait measurements are from museum skins, while others are gaps that have been imputed / estimated. This 159 includes
a handful of extinct species no longer included in the study datasets, but does not include
the four extinct species that are included in AVONET (i.e. the previous file).

• CONS_TREE – the consensus phylogeny, with 10,154 tips. Consensus tree built using multiple Jetz et al. (2012) BirdTree trees.




