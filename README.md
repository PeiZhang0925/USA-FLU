# USA-FLU
   
This repository contains all the codes and data used to produce the results of the USA-FLU paper.
  
## Data   
  
The Data folder contains all the data used in this study.  
  
(1). **usa_Flu_P_proxy_Data_B.rda**: this data set contains information of temperature, absolute humidity, ozone, and flu in flu season (October ~ May) from 2010 to 2015. This is the main data used in analyses.    
(2). **usa_Flu_P_proxy_Data_full_B.rda**: this is the data used to produce surrogate data of environmental variables in MFI significant test part. It includes temperature, absolute humidity, and ozone data from January to December from 2010 to 2015.
  
## Analyses  
  
The Analyses folder contains all the codes used in this study.  
  
(1) R_functions.R   
(2) GAM.R       
(3) PCMCI_data_processing.R   
(4) PCMCI+.py  
(5) EDM.R   
  
### Functions    
  
**R_functions.R** file is created to make the codes files more readable and userfriendly. All R functions used in GAM, PCMCI+, and EDM have been included in this file. To load these functions, just use code ```source("R_functions.R")```.
   
### 1. GAM   
  
**GAM.R** file is used for the GAM analysis. In this file, we first build GAM to obtain beta results. Then we transform beta results to risk ratio (RR) of each IQR change to be comparable to existing evidence.  
  
### 2. PCMCI+   
  
**PCMCI_data_processing.R** file is for the data processing of PCMCI+ analysis (Run this file first before PCMCI+ analysis). After running this file, we can get usaflu_raw.csv for state-specific PCMCI+ analysis and usaflu_norm.csv for overall PCMCI+ analysis.  
  
**PCMCI+.py** is used for PCMCI+ analysis. Since there is no tigramite package (the package for PCMCI+) for R, this part of analysis is done on python. In this file, we conduct overall PCMCI+ analysis for the combined data of all the states first. Then we do the same analysis for each state.    
  
### 3. EDM   
   
**EDM.R** contains all the codes for EDM analysis. The main steps in this file are listed below (also clearly shown in the R file):   
  
(1) **Calculate noise factors for surrogate data of environmental variables**. By this step, state-wise best noise factor values can be obtained for the production of surrogate data in the fifth step MFI significance test.     
(2) **Determine optimal E for the system**. E is the demension used to reconstruct the dynamic system. The E leading to the best predictive performance is selected in this step for subsequent analyses.  
(3) **Determine optimal theta for S-map**. θ is the parameter controlling nonlinearity. The greater the θ, the larger the nonlinearity. Here θ peaks the predictive performance is selected for next analyses.    
(4) **Calculate state-specific forecast improvement**. Forecast improvement can be used to evaluate the causal relationship between variables. By this step, we are going to compare the effects of temperature, absolute humidity, and ozone on flu.    
(5) **MFI significance test**. We test if the relationship between environmental variables and flu are significant in this step by compare the predictive performace of real data and surrogate data. This step takes about 2h to run.           
(6) **Effect strength**. The partial derivate of smap represents the effect stength. After testing the significance of relationships, we calculate corresponding effect strength in this step.    
  
