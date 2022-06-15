
rm(list = ls())

################################## Load packages ###################################

options(width=80)
options(stringsAsFactors=FALSE)
options(scipen = 6) # bias against scientific notation

# packages
packages=c('tidyverse','knitr','lubridate','rEDM','metap',
             'doParallel','foreach','imputeTS')
lapply(packages, require, character.only=T)

# parallel computing parameters
cores_all=detectCores()
cores=ifelse(cores_all<9,4,cores_all-2)
core_type='PSOCK'

################################## Load flu-season data and functions ###################################

## load data and functions
load("usa_Flu_P_proxy_Data_B.rda")
source("R_functions_usaflu.R")

## logit transformation of fluP
dfA=usa_Flu_P_proxy_Data_B %>% filter(state=="ak")
dfA$logitfluP <- logitTransform(dfA$fluP)
dfA$fluP <- dfA$logitfluP

num_sample=100
num_surr=1000

## normalization of variables
df_smapc=dfA %>% group_by(state) %>%
  mutate_at(3:ncol(.),normFunc) %>% ungroup()

################################# Calculate surrogate data noise factor of independent variables (all year data) #######################

load("usa_Flu_P_proxy_Data_full_B.rda")

# normalization of variables
dfB=usa_Flu_P_proxy_Data_full_B %>% filter(state=="ak")
dfB_smapc=dfB %>% group_by(state) %>%
  mutate_at(3:ncol(.),normFunc) %>% ungroup()

# select independent variables
df_flu=dfB_smapc %>% select(date,state,ah,o3,temp) %>%
  group_by(state) %>% gather(key,value,-c(date,state))

# calculate spar values for independent variables
plist=list(states=unique(df_flu$state),
           keys=c('o3','ah','temp')) %>% cross_df()

cl <- makeCluster(cores[1],type = core_type)
registerDoParallel(cl)
flu_spar=plist %>% pmap_df(fn_state_spar)
stopCluster(cl)

# calculate sd values for independent variables (additive noise factor to produce surrogate data)
plist=list(states=unique(df_flu$state),
           plts=c('o3','ah','temp')) %>% cross_df()

sd_data=plist %>% pmap_df(fn_anomaly_PNAS)

################################### Determine optimal E for the system by each state ###########################################

plist=list(data=list(df_smapc),y=c('fluP'),ST=unique(dfA$state)) %>%
  cross_df()

E_smapc_out=plist %>% pmap_df(fn_E_smapc)

E_smapc =E_smapc_out %>% filter(E %in% 2:6) %>%
  group_by(dis,ST) %>%
  filter(rho==max(rho))  %>%
  as.data.frame() 

E_smapc

for (st in unique(dfA$state)) {

  plot_E <- ggplot(E_smapc_out %>% filter(ST== st)) +
    geom_line(aes(E, rho)) +
    scale_x_continuous(1:20) +
    ggtitle(as.character(st))

  print(plot_E)
}

################################### Determine optimal theta for S-map by each state ###########################################

plist=list(data=list(df_smapc),ST=unique(dfA$state),
           dis=c('fluP'),
           theta=c(0.01, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4,
                   5, 6, 7, 8, 9)) %>% cross_df()

cl <- makeCluster(cores[1],type = core_type)
registerDoParallel(cl)
theta_out=foreach(i = 1:nrow(plist),
                  .packages = c("rEDM","tidyverse",'lubridate'),
                  .combine=rbind,
                  .inorder=FALSE)  %dopar% {
                    theta_out=plist[i,] %>% pmap_df(fn_theta_justY)
                  }
stopCluster(cl)

best_theta=theta_out %>%
  group_by(state) %>%
  filter(rho==max(rho)) %>%
  as.data.frame()

best_theta

for (st in unique(dfA$state)) {
  plot_theta <- ggplot(theta_out %>% filter(state== st)) +
    geom_line(aes(theta, rho)) +
    ggtitle(as.character(st))
  print(plot_theta)
}

################################### Calculate state-specific forecast improvement  ###########################################

plist_1=list(data=list(df_smapc),
             ST=unique(df_smapc$state),
             tp_value=-2:0,
             dis=c('fluP'),
             plt_1=c('o3','temp','ah')) %>% cross_df()
mfi_out_1 <- plist_1 %>% pmap_df(fn_mfi_1_smap)

plist_2=list(data=list(df_smapc),
             ST=unique(df_smapc$state),
             tp_value=-2:0,
             dis=c('fluP'),
             plt_1=c('o3','ah'),
             plt_2=c('temp','ah')) %>% cross_df() %>% filter(plt_1 != plt_2)
mfi_out_2 <- plist_2 %>% pmap_df(fn_mfi_2_smap)


df_0=best_theta %>% filter(dis=='fluP') %>%
  mutate(rho_0=rho) %>%
  select(state,rho_0)

df_1=mfi_out_1 %>% filter(dis=='fluP') %>%
  mutate(rho_1=round(rho,3)) %>%
  select(state,plt_1,tp_value,rho_1) %>%
  spread(plt_1,rho_1) %>%
  rename(O3='o3',AH='ah',T='temp')

df_2=mfi_out_2 %>% filter(dis=='fluP') %>%
  mutate(rho_2=round(rho,3),
         plt=paste0(substring(plt_1,1,1),substring(plt_2,1,1))) %>%
  select(state,plt,tp_value,rho_2) %>%
  spread(plt,rho_2)

df_123=left_join(df_1,df_2,by=c('state','tp_value')) %>%
  gather(plt,rho,-state,-tp_value) %>%
  left_join(df_0,by='state') %>%
  mutate(delta_rho=rho-rho_0)

df_123_fluP <- df_123

df_A2=df_123_fluP %>%
  group_by(tp_value,plt) %>%
  summarise(mean=mean(delta_rho),median=median(delta_rho),
            p_gt_0=one_sample_mean(delta_rho)$p.value)

df_A2 %>% arrange(tp_value,plt) %>%
  print(n=Inf)
################################################## MFI Significance Tets #################################################

plist=list(data=list(df_smapc),
           ST=unique(df_smapc$state),
           dis=c('fluP'),
           plt_1=c('o3','ah','temp'),
           tp_value=-2:0) %>% cross_df()

cl <- makeCluster(cores[1],type = core_type)
registerDoParallel(cl)
mfi_surr_out=plist %>% pmap_df(fn_mfi_1_smap_surr_season)

dfg2 = best_theta %>%
  select(rho_uni=rho, ST=state) %>%
  right_join(mfi_surr_out, by= 'ST') %>%
  group_by(ST, plt, tp_value)%>%
  mutate(E= 10) %>%
  mutate(grp=ifelse(i==1,'raw','surr'),
         delta_rho=rho-rho_uni)

p_function(x='rho',pH=0.99,pL=0.01)


################################################## S-map: effect strength #################################################

plist=list(data=list(df_smapc),
           ST=unique(df_smapc$state),
           dis=c('fluP'),
           plt=c('o3','ah','temp'),
           tp_value=-2:0) %>% cross_df()


cl <- makeCluster(cores[1],type = core_type)
registerDoParallel(cl)
C_out=foreach(i = 1:nrow(plist),
              .packages = c("rEDM","tidyverse","lubridate"),
              .combine=rbind,
              .export='best_theta',
              .inorder=FALSE)  %dopar% {
                C_out=plist[i,] %>% pmap_df(fn_smapc)
              }

stopCluster(cl)

tempA <- C_out %>%   group_by(tp_value,E,dis,plt) %>%
  summarise(effect=mean(effect,na.rm=T))

ggplot(tempA,aes(tp_value,effect))+geom_line()+geom_hline(yintercept = 0) +
  ylab("Effect size")+
  xlab("lag")+
  facet_grid(~plt,scales='free')+
  theme_bw() +
  theme(axis.text.x=element_text(size = 12, color= "black"),
        axis.text.y=element_text(size = 12, color= "black"),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        axis.ticks.length=unit(0.2,'cm'),
        legend.position="none",
        panel.background = element_blank(),
        plot.title = element_text(size = 18, hjust=0.5))



