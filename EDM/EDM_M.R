########################################################
## "Ozone as another environmental driver of influenza" 
##            Results Reproduce Code File
##                  
##             BY Fang
## Updates: influenza activity proxy only
## Updates: ozone data supplemented
## Updates: states of 90% complete data included followed by interpolation.
## Updates: Hawaii (HI) and Alaska (AK) included
## Updates: E selection range 1:10; pred=1:m
## Updates: E_smapc for CCM analysis
## Updates: LOOGV in E_smapc & theta selection, choose the initial peak of E and max of theta
## Updates: alpha=alpha for CCM and MFI surrogate test
## Updates: tp_value=-2:0 in CCM & MFI
## Updates: correct normalization as state-wise
## Updates: correct `one_sample_mean` function, not ignoring `tp_value` info
## Updates: state-specific MFI single delta_rho >0 surrogate test
## Updates: s-map for mfi & effect size estimation
## Updates: excluding first&last E rows in 'C_out'
##
## updates: use data version M, 26 states, WV length 236, others 261
########################################################


################################## Load packages ###################################

# opts_chunk$set(warning=FALSE,message=FALSE)
options(warn = -1)
options(width=80) 
options(stringsAsFactors=FALSE)
options(scipen = 6) # bias against scientific notation
options(digits = 3) # show fewer decimal places

# packages
packages=c('tidyverse','knitr','lubridate','rEDM','metap',
           'doParallel','foreach')

lapply(packages, require, character.only=T)


cores_all=detectCores()
cores=ifelse(cores_all<9,4,cores_all-2)
core_type='PSOCK'


################################## Load data ###################################

load('usa_Flu_M_proxy_Data.rda')

dfA=usa_Flu_M_proxy_Data 

table(dfA$state)

normFunc=function(x){(x-mean(x, na.rm = T))/sd(x, na.rm = T)}

df_smapc=dfA %>% group_by(state) %>% 
  mutate_at(3:ncol(.),normFunc) %>% ungroup()


################################# Determine optimal E for the system #################################

fn_E_smapc=function(data,ST,y){
  
  df=data %>%
    filter(state==ST) %>%
    select(date,y) %>%
    na.omit()
  
  M <- nrow(df)
  
  # lib <- c(1, floor(2/3 * M))
  # pred <- c(floor(2/3 * M) + 1, M)
  
  lib <- c(1, M)
  pred <- c(1, M)
  
  E=EmbedDimension(dataFrame=df,
                   lib=lib,
                   pred=pred,
                   columns=y, target=y,
                   maxE=20,
                   showPlot=F)
  temp=data.frame(dis=y,E,ST)
  temp}

plist=list(data=list(dfA),y=c('fluP'),ST=unique(dfA$state)) %>%
  cross_df()

# E_smapc_out=plist %>% pmap_df(fn_E_smapc)
# 
# save(E_smapc_out, file = 'E_smapc_out_M.rda')

load('E_smapc_out_M.rda')

find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pk <- unlist(pks)
  pks <- ifelse(is.null(pk)==T, which.max(x), pk)
  pks
}

E_smapc =E_smapc_out %>% filter(E %in% 1:10) %>%
  group_by(dis,ST) %>%
  mutate(row_number=1:n()) %>%
  filter(rho==max(rho))  %>% 
  as.data.frame() %>%
  select(-row_number)

E_smapc

for (st in unique(dfA$state)) {
  
  plot_E <- ggplot(E_smapc_out %>% filter(ST== st)) +
    geom_line(aes(E, rho)) +
    scale_x_continuous(1:20) +
    ggtitle(as.character(st))
  
  print(plot_E)
}

################################# Calculate parameter alpha for surrogate data #################################

load('sd_PNAS_fluP_sparAsIs_M.rda')
sd_data=sd_PNAS_fluP_sparAsIs_M

#################################### 1. CCM_OUT: perform CCM for real and surrogate data ###################################

##### 1) Perform CCM for real and surrogate data with varing alpha in surrogates generation #####

num_surr=500
num_sample=100

fn_season_ccm=function(data,ST,x,y,tp_value){
  df=data %>%
    filter(state==ST) %>%
    select(date,y,x)
  
  m=nrow(df) %>% as.data.frame()
  
  E=E_smapc[E_smapc[,'ST']==ST & E_smapc[,'dis']==y,'E']
  
  alpha=sd_data %>% filter(states==ST & plt==x) %>%
    select(sd_PNAS) %>% pull()
  
  libSize =c(m-E)
  
  set.seed(2019)
  surr_data=
    SurrogateData(unlist(df[,x]), method = "seasonal",
                  T_period = 52.18,
                  num_surr = num_surr, 
                  alpha=alpha)
  
  all_data = as.data.frame(cbind(df,surr_data))
  
  names(all_data) = c("date", y, 'T1',paste0("T", 2:(num_surr+1)))
  
  rho_surr <- data.frame(i= numeric(num_surr+1),
                         rho=numeric(num_surr+1),
                         ST=ST,dis=y,plt=x,E=E,
                         tp_value=tp_value,
                         alpha=alpha)
  
  for (i in 1:(num_surr+1)) {
    targetCol = paste("T", i, sep = "")
    ccm_out = CCM(dataFrame = all_data, E = E, Tp = tp_value,
                  columns = y,
                  target = targetCol,
                  libSizes = libSize,
                  random=F,
                  sample = num_sample,
                  seed=2019)
    col = paste(y, ":", targetCol, sep = "")
    rho_surr$rho[i] = ccm_out[1, col]
    rho_surr$i = 1:(num_surr+1)
  }
  rho_surr
}

plist=list(data=list(dfA),
           ST=unique(dfA$state),
           y=c('fluP'),
           x=c('o3',"temp","ah"),
           tp_value=-2:0) %>% cross_df()

# cl <- makeCluster(cores[1],type = core_type)
# registerDoParallel(cl)
# ccm_out=foreach(j = 1:nrow(plist),
#                 .packages = c("rEDM","tidyverse"),
#                 .combine=rbind,
#                 # .export='num_sample',
#                 .inorder=FALSE)  %dopar% {
#                   ccm_out=plist[j,] %>% pmap_df(fn_season_ccm)
#                 }
# stopCluster(cl)
# 
# save(ccm_out,file = 'ccm_out_M_max.rda')

load('ccm_out_M_max.rda')

ccm_p=ccm_out %>% 
  filter(!ST %in% c('hi',"ak")) %>%
  group_by(dis,plt,E,ST,tp_value) %>%
  summarise(p=1-ecdf(rho[i != 1])(rho[i == 1]))

fn_metap=function(var,tp_values,plts,diss){
  df=filter(var,tp_value==tp_values & plt==plts & dis==diss)
  out=allmetap(df$p, method = "sumlog") %>% as.data.frame()
  mutate(out,tp_value=tp_values,plt=plts,dis=diss)
}

plist=list(var=list(ccm_p),
           tp_values=-2:0,
           plts=c("o3","temp","ah"),
           diss=c('fluP')) %>% cross_df()

meta_ccm_p_out=plist %>% pmap_df(fn_metap)

meta_ccm_p_out



################################# 2. MFI: Multivariate Forecast improvement ##########################################

##### 1) Determine optimal θ for S-map by each state #####

fn_theta_justY=function(data, ST, dis, theta){
  
  E=E_smapc[E_smapc[,'ST']==ST & E_smapc[,'dis']==dis,'E']
  
  df=data %>%
    filter(state==ST) %>%
    #rename(date='start_d') %>%
    #mutate(date=as.Date(date,'%Y-%m-%d')) %>% 
    select(date,dis)
  
  M <- nrow(df)
  
  lib=c(1,M)
  pred=c(1,M)
  
  rho_theta = PredictNonlinear(dataFrame = df,
                               embedded = FALSE,
                               columns = dis,
                               target = dis,
                               Tp=1,
                               theta=theta,
                               lib=lib,
                               pred=pred,
                               showPlot = FALSE,
                               E = E)
  
  best_theta_df=rho_theta %>%
    mutate(dis=dis, state=ST, theta=theta)
} 

plist=list(data=list(dfA),ST=unique(dfA$state),
           dis=c('fluP'),
           theta=c(0.01, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4,
                   5, 6, 7, 8, 9)) %>% cross_df()

# cl <- makeCluster(cores[1],type = core_type)
# registerDoParallel(cl)
# 
# theta_out=foreach(i = 1:nrow(plist),
#                   .packages = c("rEDM","tidyverse"),
#                   .combine=rbind,
#                   .inorder=FALSE)  %dopar% {
#                     theta_out=plist[i,] %>% pmap_df(fn_theta_justY)
#                   }
# stopCluster(cl)
# 
# save(theta_out, file = 'theta_out_M_max.rda')

load('theta_out_M_max.rda')

best_theta=theta_out %>%
  group_by(state) %>%
  filter(rho==max(rho)) %>%
  # mutate(row_number=1:n()) %>%
  # filter(row_number==find_peaks(rho,m=0)[1]) %>%
  as.data.frame()

E_smapc1 <-E_smapc
names(E_smapc1) <- c("dis","E","rho","state" )
best_theta1 <- merge(best_theta, E_smapc1[,c("E","state")], by="state")

best_theta1

for (st in unique(dfA$state)) {
  
  plot_theta <- ggplot(theta_out %>% filter(state== st)) +
    geom_line(aes(theta, rho)) +
    # scale_x_continuous(1:20) +
    ggtitle(as.character(st))
  
  print(plot_theta)
}


##### 3) Calculate state-specific forecast improvement #####

fn_mfi_0_smap=function(data,ST,dis){
  E=E_smapc[E_smapc$ST==ST & E_smapc$dis==dis,"E"]
  df=data %>% 
    filter(state==ST) %>%
    select(date,dis) %>%
    mutate(date=as.Date(date,'%Y-%m-%d')) 
  
  M <- nrow(df)
  
  embed_1=Embed(dataFrame = df, E = E, tau = -1, columns = dis )
  
  
  dataFrame_0 = cbind(df[E:M, 'date'],df[E:M, dis],
                      embed_1[E:M, 1:E]) %>% as.data.frame() 
  
  names(dataFrame_0)=c('date',dis,letters[1:E])
  columns_0 = paste(letters[1:E],sep =' ') 
  
  m = nrow(dataFrame_0)
  
  theta=best_theta[best_theta$state==ST & best_theta$dis==dis,"theta"]
  
  smap_0 = SMap(dataFrame = dataFrame_0,
                embedded = TRUE,
                columns = columns_0,
                target = dis,
                lib = c(1, m),
                pred = c(1, m),
                theta=theta,
                # knn = E+4,
                Tp = 1)
  
  out_0=compute_stats(smap_0$predictions$Observations,
                      smap_0$predictions$Predictions ) 
  
  out=rbind(out_0) %>%
    mutate(state=ST,dis=dis)
  out
}

# data=df_smapc; ST='ak'; dis='logitfluP'; plt_1='o3'; plt_2='ah'; tp_value=0

fn_mfi_1_smap=function(data,ST,plt_1,dis,tp_value){
  E=E_smapc[E_smapc$ST==ST & E_smapc$dis==dis,"E"]
  df=data %>% 
    filter(state==ST) %>%
    # rename(date='start_d') %>%
    mutate(plt_1=lag(.data[[plt_1]],-tp_value)) %>% 
    filter(!(is.na(plt_1))) %>%
    select(date,dis,plt_1) %>%
    mutate(date=as.Date(date,'%Y-%m-%d')) 
  
  M <- nrow(df)
  
  embed_1=Embed(dataFrame = df, E = E, tau = -1, columns = dis )
  
  dataFrame_1 = cbind(df[E:M, 'date'],df[E:M, dis],
                      embed_1[E:M, 1:E],df[E:M, 'plt_1']  ) %>% 
    as.data.frame()
  
  names(dataFrame_1)=c('date',dis,letters[1:E],plt_1)
  columns_1 = paste(paste(letters[1:E],collapse =' '), plt_1 ,
                    sep=' ')
  
  m = nrow(dataFrame_1)
  
  theta=best_theta[best_theta$state==ST & best_theta$dis==dis,"theta"]
  
  smap_1 = SMap(dataFrame = dataFrame_1,
                embedded = TRUE,
                columns = columns_1,
                target = dis,
                lib = c(1, m),
                pred = c(1, m),
                theta=theta,
                # knn = E+4,
                Tp = 1)
  
  out_1=compute_stats(smap_1$predictions$Observations,
                      smap_1$predictions$Predictions )  %>%
    mutate(mfi_012=1)
  
  out=rbind(out_1) %>%
    mutate(state=ST,plt_1=plt_1,dis=dis,
           tp_value=tp_value)
  out
}

fn_mfi_2_smap=function(data,ST,plt_1,plt_2,dis,tp_value){
  E=E_smapc[E_smapc$ST==ST & E_smapc$dis==dis,"E"]
  df=data %>% 
    filter(state==ST) %>%
    # rename(date='start_d') %>%
    mutate(plt_1=lag(.data[[plt_1]],-tp_value)) %>% 
    mutate(plt_2=lag(.data[[plt_2]],-tp_value)) %>% 
    filter(!(is.na(plt_1))) %>%
    filter(!(is.na(plt_2))) %>%
    select(date,dis,plt_1,plt_2) %>%
    mutate(date=as.Date(date,'%Y-%m-%d')) 
  
  M <- nrow(df)
  
  embed_1=Embed(dataFrame = df, E = E, tau = -1, columns = dis )
  
  dataFrame_2 = cbind(df[E:M, 'date'],df[E:M, dis],
                      embed_1[E:M, 1:E],df[E:M, 'plt_1'],
                      df[E:M, 'plt_2']) %>% as.data.frame()
  
  names(dataFrame_2)=c('date',dis,letters[1:E],plt_1,plt_2)
  columns_2 = paste(paste(letters[1:E],collapse =' '), plt_1,plt_2 ,
                    sep=' ')
  
  m = nrow(dataFrame_2)
  
  theta=best_theta[best_theta$state==ST & best_theta$dis==dis,"theta"]
  
  smap_2 = SMap(dataFrame = dataFrame_2,
                embedded = TRUE,
                columns = columns_2,
                target = dis,
                lib = c(1, m),
                pred = c(1, m),
                theta=theta,
                # knn = E+2,
                Tp = 1)
  
  out_2=compute_stats(smap_2$predictions$Observations,
                      smap_2$predictions$Predictions )  
  
  out=rbind(out_2) %>%
    mutate(state=ST,plt_1=plt_1,plt_2=plt_2,dis=dis,
           tp_value=tp_value)
  out
}


plist_0=list(data=list(df_smapc),
             ST=unique(df_smapc$state),
             dis=c('fluP')) %>%
  cross_df()
mfi_out_0 <- plist_0 %>% pmap_df(fn_mfi_0_smap)

plist_1=list(data=list(df_smapc),
             ST=unique(df_smapc$state),
             tp_value=-2:0,
             dis=c('fluP'),
             plt_1=c('o3','temp','ah')) %>%
  cross_df()
mfi_out_1 <- plist_1 %>% pmap_df(fn_mfi_1_smap)

plist_2=list(data=list(df_smapc),
             ST=unique(df_smapc$state),
             tp_value=-2:0,
             dis=c('fluP'),
             plt_1=c('o3','ah'),
             plt_2=c('temp','ah')) %>%
  cross_df() %>% filter(plt_1 != plt_2)
mfi_out_2 <- plist_2 %>% pmap_df(fn_mfi_2_smap)


df_0=mfi_out_0 %>% filter(dis=='fluP') %>%
  mutate(rho_0=round(rho,3)) %>%
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
  #left_join(df_3,by=c('state','tp_value')) %>%
  gather(plt,rho,-state,-tp_value) %>%
  left_join(df_0,by='state') %>%
  mutate(delta_rho=rho-rho_0)

df_123_fluP <- df_123

# save(df_123_fluP, file = 'mfi_out_M_max.rda')

load('mfi_out_M_max.rda')

################################################## MFI Significance Tets #################################################

num_surr=500
num_sample=100

plist=list(data=list(df_smapc),
           ST=unique(df_smapc$state),
           dis=c('fluP'),
           plt_1=c('o3','ah','temp'),
           tp_value=-2:0) %>% cross_df()

##### 1) Perform surrogate test with varying alpha in surrogates generation #####

cl <- makeCluster(cores[1],type = core_type)
registerDoParallel(cl)
fn_mfi_1_smap_surr_season <-
  function(data,ST,plt_1,dis,tp_value){
    
    E=E_smapc[E_smapc$ST==ST & E_smapc$dis==dis,"E"]
    df <- data %>% filter(state==ST) %>%
      mutate(plt=lag(.data[[plt_1]],-tp_value)) %>%
      filter(!(is.na(plt))) %>% select(date,dis,plt) %>%
      mutate(date=as.Date(date,'%Y-%m-%d'))%>%
      na.omit()
    
    alpha=sd_data %>% filter(states==ST & plt==plt_1) %>%
      select(sd_PNAS) %>% pull()
    
    set.seed(2019)
    surr_data <- SurrogateData(ts(df[,'plt']),
                               method = "seasonal",
                               T_period = 52.18,
                               num_surr = num_surr,
                               alpha = alpha)
    
    all_data <- as.data.frame(cbind(df,surr_data))
    names(all_data) <- c("date", dis, 'T1',paste0("T", 2:(num_surr+1)))
    M <- nrow(all_data)
    out <- foreach(i = 1:(num_surr+1),
                   .packages = c("rEDM","tidyverse"),
                   .export=c('num_sample','num_surr','best_theta'),
                   .combine=rbind,
                   .inorder=FALSE)  %dopar% {
                     targetCol <- paste("T", i, sep = "")
                     embed_1 <- Embed(dataFrame = all_data, E = E, tau = -1, columns = dis )
                     dataFrame_1 <- cbind(all_data[E:M, 'date'],all_data[E:M, dis],embed_1[E:M, 1:E],
                                          all_data[E:M, targetCol]) %>% as.data.frame()
                     names(dataFrame_1) <- c('date',dis,letters[1:E],plt_1)
                     columns_1 <- paste(paste(letters[1:E],collapse =' '), plt_1 , sep=' ')
                     theta=best_theta[best_theta$state==ST & best_theta$dis==dis,"theta"]
                     smap_1 <- SMap(dataFrame = dataFrame_1,embedded = TRUE,
                                    columns = columns_1,target = dis,
                                    lib = c(1, nrow(dataFrame_1)),
                                    pred = c(1, nrow(dataFrame_1)),
                                    theta = theta,Tp = 1)
                     out_1=compute_stats(smap_1$predictions$Observations,
                                         smap_1$predictions$Predictions) %>% mutate(mfi_012=1)
                     out=rbind(out_1) %>% mutate(i=i,ST=ST,plt=plt_1,dis=dis,tp_value=tp_value, 
                                                 alpha=alpha)
                     out
                   }
    out
  }



# mfi_surr_out=plist %>% pmap_df(fn_mfi_1_smap_surr_season)
#  
# save(mfi_surr_out,file = 'mfi_surr_out_M_max.rda')

load('mfi_surr_out_M_max.rda')

# mfi_surr_out %>%
#   filter(i==1) %>%
#   select(state=ST, tp_value, plt, rho) %>%
#   mutate(rho=round(rho,3))%>%
#   spread(plt,rho)%>%
#   left_join(df_1, by=c('state', 'tp_value'))

dfg2 = best_theta %>%
  select(rho_uni=rho, ST=state) %>%
  right_join(mfi_surr_out, by= 'ST') %>%
  group_by(ST, plt, tp_value)%>%
  mutate(E= 10) %>%
  mutate(grp=ifelse(i==1,'raw','surr'),
         delta_rho=rho-rho_uni) 

dfg2.p=dfg2 %>%
  group_by(dis,plt,E,ST,tp_value) %>%
  summarise(p=1-ecdf(delta_rho[i != 1])(delta_rho[i == 1])) %>%
  arrange(dis,plt)

print(dfg2.p, n=Inf)  # comparing raw and surrogate time series

####################################
####################################
fn_metap=function(var,tp_values,plts,diss){
  df=filter(var,tp_value==tp_values & plt==plts & dis==diss)
  out=allmetap(df$p, method = "sumlog") %>% as.data.frame()
  mutate(out,tp_value=tp_values,plt=plts,dis=diss)
}

plist=list(var=list(dfg2.p),
           tp_values=-2:0,
           plts=c("o3","temp","ah"),
           diss=c('fluP')) %>% cross_df()

meta_mfi_p_out=plist %>% pmap_df(fn_metap)

meta_mfi_p_out

plist25=list(var=list(dfg2.p %>% filter(!ST %in% c("ak","hi"))),
             tp_values=-2:0,
             plts=c("o3","temp","ah"),
             diss=c('fluP')) %>% cross_df()

meta_mfi_p_out25=plist25 %>% pmap_df(fn_metap)

meta_mfi_p_out25
##############################  3. S-map: effect strength ############################

# data=df_smapc; ST='ga'; dis='fluP'; plt='o3'; tp_value=-2; theta=0.1

fn_smapc=function(data,ST,plt,dis,tp_value){
  
  E=E_smapc[E_smapc$ST==ST & E_smapc$dis==dis,"E"]
  
  df=data %>% filter(state==ST) %>%
    mutate(plt_tp=lag(.data[[plt]],-tp_value)) %>%
    filter(!(is.na(plt_tp))) %>%
    select(date,all_of(dis),plt_tp) %>%
    na.omit()
  
  M <- nrow(df)
  
  embed_1=Embed(dataFrame = df, E = E, tau = -1, columns = dis )
  
  dataFrame = cbind(df[E:M, 'date'],df[E:M, dis],
                    embed_1[E:M, 1:(E-1)], df[E:M, 'plt_tp']  ) %>%
    as.data.frame()
  
  names(dataFrame)=c('date',dis,letters[1:(E-1)],plt)
  
  columns = paste(paste(letters[1:(E-1)],collapse =' '), plt ,
                  sep=' ')
  
  m <- nrow(dataFrame)
  
  lib=c(1,m)
  pred=c(E + 1, m - E)
  
  theta = best_theta[best_theta$state==ST & best_theta$dis==dis, "theta"]
  
  smap = SMap(dataFrame = dataFrame,
              embedded = TRUE,
              columns = columns,
              target = dis,
              lib = lib,
              pred=pred,
              theta=theta,
              Tp = 1,
              E = E)
  smapc_df=smap$coefficients[c(1,2+E)]
  
  names(smapc_df)=c('date','effect')
  smapc_df=smapc_df %>%
    mutate(date=lubridate::as_date(date,
                                   origin = lubridate::origin)) %>%
    mutate(dis=dis,ST=ST, plt=plt,E=E,tp_value=tp_value)
}

# fn_smapc(data = df_smapc, ST='ga',dis=c('fluP'),
#         plt=c('o3'), tp_value = 0)

plist=list(data=list(df_smapc),
           ST=unique(df_smapc$state),
           dis=c('fluP'),
           plt=c('o3','ah','temp'),
           tp_value=-2:0) %>% cross_df()


# cl <- makeCluster(cores[1],type = core_type)
# registerDoParallel(cl)
# C_out=foreach(i = 1:nrow(plist),
#               .packages = c("rEDM","tidyverse","lubridate"),
#               .combine=rbind,
#               .export='best_theta',
#               .inorder=FALSE)  %dopar% {
#                 C_out=plist[i,] %>% pmap_df(fn_smapc)
#               }
# 
# stopCluster(cl)
# 
# save(C_out, file = 'C_out_M_max.rda')

load('C_out_M_max.rda')

head(C_out %>% arrange(plt),n=10)
C_out %>% group_by(tp_value,ST,dis,plt) %>%
  slice(E+1:E+2)%>% arrange(plt,ST,tp_value)

tempA=C_out %>% group_by(tp_value,ST,E,dis,plt) %>%
  # slice((E+1):(nrow(.)-E)) %>%
  summarise(effect=mean(effect,na.rm=T))

tempB = tempA  %>% ungroup() %>%
  select(-E,-dis)

mfip_effect_each = tempB %>% 
  left_join(dfg2.p %>% ungroup() %>% select(ST, plt,tp_value, p), by=c('ST', 'plt','tp_value')) %>%
  arrange(ST,plt,tp_value)

tempC = tempB %>%
  group_by(plt, tp_value) %>%
  dplyr::summarise(n = n(),
                   effect =mean(effect,na.rm=T)) %>%
  print(n= 100)


######################## 3 Methods Table ################################
# EDM_each <- mfip_effect_each %>% mutate(EDM.effect=effect,mfi.p=p) %>%
#   select(ST,plt,tp_value,EDM.effect,mfi.p)
# 
# EDM_meta <- meta_mfi_p_out %>%
#   mutate(ST='Overall') %>%
#   rename(mfi.p = p) %>%
#   right_join(tempC, by=c('plt','tp_value')) %>%
#   rename(EDM.effect= effect)%>%
#   select(ST, plt, tp_value, EDM.effect, mfi.p)
# 
# 
# 
# EDM <- rbind(EDM_each,EDM_meta) %>% mutate(EDM.sig=ifelse(mfi.p<=0.05,1,0)) 
# EDM$tp_value <- abs(EDM$tp_value)
# 
# # 
# # EDM_B <- EDM_A %>%
# #   rename(mfi.p = p) %>%
# #   left_join(ccm_p, by=c('plt','ST','tp_value')) %>%
# #   rename(ccm.p = p, EDM.effect= effect) %>%
# #   select(ST, plt, tp_value, EDM.effect, mfi.p, ccm.p)%>%
# #   rbind(meta_edm) %>%
# #   mutate(ST=case_when(ST == "al" ~ "Alabama", ST == "ar" ~ "Arkansas",
# #                       ST == "az" ~ "Arizona", ST == "ca" ~ "California",
# #                       ST == "co" ~ "Colorado", ST == "ct" ~ "Connecticut",
# #                       ST == "de" ~ "Delaware", ST == "ga" ~ "Georgia",
# #                       ST == "il" ~ "Illinois", ST == "ky" ~ "Kentucky",
# #                       ST == "la" ~ "Louisiana", ST == "md" ~ "Maryland",
# #                       ST == "mn" ~ "Minnesota", ST == "mo" ~ "Missouri",
# #                       ST == "mt" ~ "Montana", ST == "ny" ~ "New York",
# #                       ST == "oh" ~ "Ohio", ST == "ok" ~ "Oklahoma",
# #                       ST == "or" ~ "Oregon", ST == "pa" ~ "Pennsylvania",
# #                       ST == "sc" ~ "South Carolina", ST == "sd" ~ "South Dakota",
# #                       ST == "tx" ~ "Texas", ST == "ut" ~ "Utah",
# #                       ST == "va" ~ "Virginia", ST == "wa" ~ "Washington",
# #                       ST == "wi" ~ "Wisconsin", ST == "wv" ~ "West Virginia",
# #                       ST == "hi" ~ "Hawaii", ST == "ak" ~ "Alaska",
# #                       TRUE ~ ST)) %>%
# #   rename(state=ST)%>%
# #   mutate(lag=case_when(tp_value==0 ~ 0,
# #                        tp_value==-1 ~ 1,
# #                        tp_value==-2 ~ 2))%>%
# #   select(state, plt, lag, EDM.effect, mfi.p, ccm.p)
# 
# ## -----
# 
# PCMCI <- read.csv('final_res_pcmci.csv')
# 
# PCMCI <- PCMCI %>%
#   rename(PCMCI.effect=effect, pcmci.sig=link,ST=state,tp_value=lag)%>%
#   select(ST,plt,tp_value,PCMCI.effect,pcmci.sig)
# PCMCI_EDM <- PCMCI %>% left_join(EDM,by=c("ST","plt","tp_value")) %>% as.data.frame()
# PCMCI_EDM$mfi.p <- unlist(PCMCI_EDM$mfi.p) %>% round(5)
# write.csv(PCMCI_EDM,file = "PCMCI_EDM_final.csv")
# PCMCI_EDM <- read.csv( "PCMCI_EDM_final.csv")
# GAM <- read.csv('GAMres.csv')
# names(GAM)
# names(PCMCI_EDM)
# GAM <- GAM %>% rename(GAM.effect=effect, ST=state,tp_value=lag) %>% 
#   select(ST,plt,tp_value,GAM.effect,GAM.sig)
# 
# PCMCI_EDM_GAM <- PCMCI_EDM %>%
#   left_join(GAM, by=c("ST","plt","tp_value")) %>% select(-mfi.p)
# PCMCI_EDM_GAM$PCMCI_EDM_GAM <- PCMCI_EDM_GAM$pcmci.sig+PCMCI_EDM_GAM$EDM.sig+PCMCI_EDM_GAM$GAM.sig
# 
# PCMCI_EDM_GAM$PCMCI_EDM <- PCMCI_EDM_GAM$pcmci.sig+PCMCI_EDM_GAM$EDM.sig
# write.csv(PCMCI_EDM_GAM, "PCMCI_EDM_GAM_final.csv")
# PCMCI_EDM_GAM$PCMCI_EDM <- as.factor(PCMCI_EDM_GAM$PCMCI_EDM)
# library(ggrepel)
# ggplot(PCMCI_EDM_GAM,aes(PCMCI.effect,EDM.effect,col=PCMCI_EDM))+geom_point()+
#   geom_text_repel(aes(PCMCI.effect,EDM.effect, label = ST),size=3) +
#   facet_grid(plt~tp_value)
# 
# 
# # save(PCMCI_EDM_GAM, file = 'PCMCI_EDM_GAM.rda')
# #
# # write.csv(PCMCI_EDM_GAM, file = 'PCMCI_EDM_GAM.csv')
# 
# # roundFunc <- function(x){round(x, 2)}
# # 
# # PCMCI_EDM_GAM_o3 <- PCMCI_EDM_GAM %>%
# #   filter(plt=='o3') %>%
# #   select(-X) %>%
# #   mutate_at(4:ncol(.), roundFunc)
# # 
# # # write.csv(PCMCI_EDM_GAM_o3, file = 'PCMCI_EDM_GAM_o3_H.csv')
# # head(PCMCI_EDM_GAM_o3, n=20)
# # 
# # GAM_beta <- read.csv('tableres_GAM_o3_B_beta.csv')
# # 
# # PCMCI_EDM_GAM_beta <- PCMCI_EDM %>%
# #   left_join(GAM_beta, by=c('state','plt','lag'))%>%
# #   mutate(mfi.p=as.numeric(mfi.p), ccm.p=as.numeric(ccm.p))%>%
# #   rename(GAM.effect=effect, gam.p=P_value)
# # 
# # PCMCI_EDM_GAM_o3_beta <- PCMCI_EDM_GAM_beta %>%
# #   filter(plt=='o3') %>%
# #   select(-X) %>%
# #   mutate_at(4:ncol(.), roundFunc)
# # 
# # # write.csv(PCMCI_EDM_GAM_o3_beta, file = 'PCMCI_EDM_GAM_o3_K_beta.csv')
# # head(PCMCI_EDM_GAM_o3_beta, n=20)
# # # 
# # # # save.image("2021flu_o3_K.rda")
# # # 
# # para_slct <- E_smapc %>% select(E, state=ST) %>%
# #   left_join(best_theta, by='state') %>%
# #   select(state, E, theta) %>%
# #   rename(ST=state) %>%
# #   mutate(ST=case_when(ST == "al" ~ "Alabama", ST == "ar" ~ "Arkansas",
# #                       ST == "az" ~ "Arizona", ST == "ca" ~ "California",
# #                       ST == "co" ~ "Colorado", ST == "ct" ~ "Connecticut",
# #                       ST == "de" ~ "Delaware", ST == "ga" ~ "Georgia",
# #                       ST == "il" ~ "Illinois", ST == "ky" ~ "Kentucky",
# #                       ST == "la" ~ "Louisiana", ST == "md" ~ "Maryland",
# #                       ST == "mn" ~ "Minnesota", ST == "mo" ~ "Missouri",
# #                       ST == "mt" ~ "Montana", ST == "ny" ~ "New York",
# #                       ST == "oh" ~ "Ohio", ST == "ok" ~ "Oklahoma",
# #                       ST == "or" ~ "Oregon", ST == "pa" ~ "Pennsylvania",
# #                       ST == "sc" ~ "South Carolina", ST == "sd" ~ "South Dakota",
# #                       ST == "tx" ~ "Texas", ST == "ut" ~ "Utah",
# #                       ST == "va" ~ "Virginia", ST == "wa" ~ "Washington",
# #                       ST == "wi" ~ "Wisconsin", ST == "wv" ~ "West Virginia",
# #                       ST == "hi" ~ "Hawaii", ST == "ak" ~ "Alaska",
# #                       TRUE ~ ST))
# # 
# # 
# # print(para_slct)
# # 
# # # write.csv(para_slct, file = 'para_slct_EDM_K.csv')
