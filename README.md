# USA-FLU
   
This repository contains all the codes and data used to produce the results of the USA-FLU paper.
  
## Data   
  
The Data folder contains all the data used in this study.  
  
1) **usa_Flu_P_proxy_Data_B.rda**: this data set contains information of temperature, absolute humidity, ozone, and flu in flu season (October ~ May) from 2010 to 2015. This is the main data used in analyses.    
2) **usa_Flu_P_proxy_Data_full_B.rda**: this is the data used to produce surrogate data of environmental variables in MFI significant test part. It includes temperature, absolute humidity, and ozone data from January to December from 2010 to 2015.
  
## Analyses  
  
1) GLM folder: Generalized Linear model, GLM_N3.R.  
2) EDM folder: Empirical dynamic modeling, EDM_N3.R (sd_PNAS_fluP_sparAsIs_N3_full.R is used to calculate spar, run it first).   
3) PCMCI+ folder: PCMCI_N3.py. data_processing_N3.R file is used to produce the data to be used in state-specific (usaflu_raw_M.csv) and nation-specific (usaflu_norm_M.csv) PCMCI+ analyses.        
