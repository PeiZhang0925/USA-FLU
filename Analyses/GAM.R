rm(list = ls())

################################## Load packages ###################################

# opts_chunk$set(warning=FALSE,message=FALSE)
options(warn = -1)
options(width=80)
options(stringsAsFactors=FALSE)
options(scipen = 6) # bias against scientific notation
options(digits = 3) # show fewer decimal places

library(tidyverse); library(lubridate); library(dlnm); library(splines);
library(tsModel); library(gnm);library(ggpubr);library(metafor);
library(mgcv); library(imputeTS)

################################## Load data and functions ###################################

load("usa_Flu_P_proxy_Data_B.rda")
source("R_functions.R")

dfA=usa_Flu_P_proxy_Data_B %>%
  mutate(time=as.numeric(date),
         year=year(date),
         month=month(date),
         week=week(date),
         ym=paste(year(date), month(date),sep = ""))

#################### GAM analysis for the three variables (get beta) ###################

plist=list(dat=list(dfA), Y='fluP',x1='o3',xlag=0:2,
           cityname=unique(dfA$state)) %>% cross_df()
out_o3  =  plist %>% pmap_df(fgam_o3)

plist=list(dat=list(dfA), Y='fluP',x1='ah',xlag=0:2,
           cityname=unique(dfA$state)) %>% cross_df()
out_ah  =  plist %>% pmap_df(fgam_ah)

plist=list(dat=list(dfA), Y='fluP',x1='temp',xlag=0:2,
           cityname=unique(dfA$state)) %>% cross_df()
out_temp  =  plist %>% pmap_df(fgam_temp)

#################### meta analysis of GAM beta results ###################

gam_o3=fmeta(out_o3)

gam_ah=fmeta(out_ah)

gam_temp=fmeta(out_temp)

#################### transform GAM beta results to effect size of each IQR change ###################

# calculate overall and state-specific IQR 
df_iqr_all=dfA %>%
  dplyr::summarize(o3=IQR(o3),
                   ah=IQR(ah),
                   temp=IQR(temp)) %>%
  pivot_longer(cols=o3:temp,names_to='plt',values_to='iqr') %>%
  mutate(state='ALL',
         iqr_all=iqr)

df_iqr_ST <- dfA %>% group_by(state) %>%
  dplyr::summarize(o3=IQR(o3),
                   ah=IQR(ah),
                   temp=IQR(temp)) %>%
  gather(plt, iqr, -state) %>%
  mutate(state=toupper(state))

df_iqr <- rbind(df_iqr_ST, df_iqr_all%>%select(state,plt,iqr))

df_iqr_plusALL <- df_iqr %>% left_join(df_iqr_all%>%select(plt,iqr_all), by='plt')

df_iqr_plusALL %>%  print(n=Inf)

# transform beta to effect size of each IQR change
df1=gam_o3  %>% as.data.frame() %>% mutate(plt='o3')
df2=gam_ah  %>% as.data.frame() %>% mutate(plt='ah')
df3=gam_temp %>% as.data.frame() %>% mutate(plt='temp')

df_gam_final=rbind(df1,df2,df3) %>%
  left_join(df_iqr_plusALL,by=c("state",'plt')) %>%
  mutate(RR= exp(beta), RRlow= exp(betalow), RRhigh=exp(betahigh)) %>%
  mutate(Size=beta*iqr,SizeL=betalow*iqr,SizeH=betahigh*iqr) %>%
  mutate(RR_iqr=exp(Size),RR_iqr_l=exp(SizeL),RR_iqr_h=exp(SizeH)) %>%
  select(state, plt, lag, beta, betalow, betahigh, RR, RRlow, RRhigh,
         Size, SizeL, SizeH, RR_iqr,RR_iqr_l,RR_iqr_h, gam.p=p, type) %>%
  mutate(sig=ifelse(gam.p<0.05, 'sig', 'non_sig'),
         sig=factor(sig, levels=c("non_sig","sig")))

df_all = df_gam_final %>% filter(state=='ALL')

df_all


