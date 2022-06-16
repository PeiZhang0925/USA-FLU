# USA-FLU
   
This repository contains all the codes and data used to produce the results of the USA-FLU paper.
  
## Data   
  
The Data folder contains all the data used in this study.  
  
**1. usa_Flu_P_proxy_Data_B.rda**: this data set contains information of temperature, absolute humidity, ozone, and flu in flu season (October ~ May) from 2010 to 2015. This is the main data used in analyses.    
**2. usa_Flu_P_proxy_Data_full_B.rda**: this is the data used to produce surrogate data of environmental variables in MFI significant test part. It includes temperature, absolute humidity, and ozone data from January to December from 2010 to 2015.
  
## Analyses  
  
The Analyses folder contains all the codes used in this study.  
  
1) R_functions.R   
2) PCMCI_data_processing.R   
3) GAM.R  
4) PCMCI+.py  
5) EDM.R   
  
### Functions    
  
**R_functions.R** file is created to make the codes files more readable and userfriendly. All R functions used in GAM, PCMCI+, and EDM have been included in this file. To load these functions, just use code ```source("R_functions.R")```.
   
### 1. GAM   
  
**GAM.R** file is used for the GAM analysis. In this file, we first build GAM to obtain beta results. Then we transform beta results to risk ratio (RR) of each IQR change to be comparable to existing evidence.  
  
### 2. PCMCI+   
  
**PCMCI_data_processing.R** file is for the data processing of PCMCI+ analysis (Run this file first before PCMCI+ analysis). After running this file, we can get usaflu_raw.csv for state-specific PCMCI+ analysis and usaflu_norm.csv for overall PCMCI+ analysis. 
  
### 3. EDM   
   
1) GLM folder: Generalized Linear model, GLM_N3.R.  
2) EDM folder: Empirical dynamic modeling, EDM_N3.R (sd_PNAS_fluP_sparAsIs_N3_full.R is used to calculate spar, run it first).   
3) PCMCI+ folder: PCMCI_N3.py. data_processing_N3.R file is used to produce the data to be used in state-specific (usaflu_raw_M.csv) and nation-specific (usaflu_norm_M.csv) PCMCI+ analyses.        
