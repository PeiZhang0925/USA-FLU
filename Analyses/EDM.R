################################## Load packages, data, and functions ###################################
# packages
packages=c('tidyverse','knitr','lubridate','rEDM','metap',
           'doParallel','foreach','imputeTS', "kableExtra")
lapply(packages, require, character.only=T)

# parallel computing parameters
cores_all=detectCores()
cores=ifelse(cores_all<9,4,cores_all-2)
core_type='PSOCK'

num_sample=100
num_surr=1000

## load data and functions
load("Data_fluseason.rda")
df <- usa_Flu_P_proxy_Data_B
# df <- read.csv("state_ky.csv")[,-1]
df$date <- as.Date(df$date )

## logit transformation of fluP
dfA=df 
logitTransform <- function(p) { log(p/(1-p)) }
dfA$logitfluP <- logitTransform(dfA$fluP)
dfA$fluP <- dfA$logitfluP

## normalization of variables
normFunc=function(x){(x-mean(x, na.rm = T))/sd(x, na.rm = T)}

df_smapc=dfA %>% mutate_at(3:ncol(.),normFunc) %>% ungroup()

################### Calculate noise factors for surrogate data of environmental variables (the whole-year data) #######################

# whole-year environmental data
load("Data_wholeyear.rda")
df <- usa_Flu_P_proxy_Data_full_B
# df <- read.csv("state_ky_full.csv")[,-1] %>% as.data.frame()
df$date <- as.Date(df$date)

# normalization of variables
dfB=df
dfB_smapc=dfB %>% group_by(state) %>%
  mutate_at(3:ncol(.),normFunc) %>% ungroup()

#select variables
df_flu=dfB_smapc %>% select(date,state,ah,o3,temp) %>%
  group_by(state) %>% gather(key,value,-c(date,state))

df_flu <- dfB_smapc %>%
  select(date, state, ah, o3, temp) %>%
  mutate_at(vars(ah, o3, temp), ~{attributes(.) <- NULL; .}) %>%
  group_by(state) %>%
  gather(key, value, -c(date, state))

# calculate spar values for independent variables
fn_state_spar=function(states,keys){
  
  splineres <- function(spar){
    res <- rep(0, length(x))
    for (i in 1:length(x)){
      mod <- smooth.spline(x[-i], y[-i], spar = spar)
      res[i] <- predict(mod, x[i])$y - y[i]
    }
    return(sum(res^2))
  }
  
  x=df_flu %>% ungroup() %>%
    filter(state==states & key==keys) %>%
    select(date) %>% pull() %>% as.numeric()
  
  y=df_flu %>% ungroup() %>%
    filter(state==states & key==keys) %>%
    select(value)  %>% pull() %>% as.numeric()
  
  spars <- seq(0, 1.5, by = 0.1)
  ss <- rep(0, length(spars))
  
  ss=foreach(i = 1:length(spars),
             # .packages = c("rEDM","tidyverse"),
             # .export=c('num_sample','num_surr'),
             .combine=rbind,
             .inorder=FALSE)  %dopar% {
               targetCol = paste("T", i, sep = "")
               ss[i] <- splineres(spars[i])
             }
  
  spar=spars[which.min(ss)]
  data.frame(state=states,plt=keys,spar)
}


plist=list(states=unique(df_flu$state),
           keys=c('o3','ah','temp')) %>% cross_df()

cl <- makeCluster(cores[1],type = core_type)
registerDoParallel(cl)
flu_spar=plist %>% pmap_df(fn_state_spar)
stopCluster(cl)

# calculate sd values (alpha) for independent variables (additive noise factor to produce surrogate data)
yearday_anom <- function(t,x,spars){
  # t: date formatted with POSIXt
  # x: time-series values to compute seasonal mean and anomaly
  doy <- as.numeric(strftime(t, format = "%j"))
  I_use <- which(!is.na(x))
  # create time indices to use for smoothing, replicating data to "wrap around"
  doy_sm <- rep(doy[I_use],3) + rep(c(-366,0,366),each=length(I_use))
  x_sm <- rep(x[I_use],3)
  xsp <- smooth.spline(doy_sm, y = x_sm, w = NULL,
                       spar = spars, cv = NA,
                       all.knots = TRUE,keep.data = TRUE, df.offset = 0)
  xbar <- data.frame(t=t,doy=doy) %>%
    left_join(data.frame(doy=xsp$x,xbar=xsp$y),by='doy') %>%
    select(xbar)
  out = data.frame(t=t,mean=xbar,anomaly=(x - xbar))
  names(out) <- c('date','mean','anomaly')
  return(out)
}

fn_anomaly_PNAS=function(states,plts){
  vec_t=df_flu %>%
    ungroup() %>%
    filter(state==states & key==plts) %>%
    select(date) %>% pull()
  
  vec_x=df_flu %>%
    ungroup() %>%
    filter(state==states & key==plts) %>%
    select(value) %>% pull()
  
  spars=flu_spar %>% filter(state==states & plt==plts) %>%
    select(spar) %>% pull()
  
  df_9=yearday_anom(vec_t,vec_x,spars)
  sd=sd(df_9$anomaly,na.rm=TRUE)
  data.frame(states,plt=plts,spar=spars,sd_PNAS=sd)
}

plist=list(states=unique(df_flu$state),
           plts=c('o3','ah','temp')) %>% cross_df()

sd_data=plist %>% pmap_df(fn_anomaly_PNAS)
sd_data

fn_surr_data=function(data,ST,plts, tp_value){
  df=data %>%
    filter(state==ST) %>%
    mutate(plt=lag(.data[[plts]],-tp_value)) %>%
    filter(!(is.na(plt))) %>%
    select(date,"plt")
  
  alpha=sd_data %>% filter(states==ST & plt==plts) %>%
    select(sd_PNAS) %>% pull()
  
  set.seed(2019)
  surr_data=
    SurrogateData(unlist(df[,"plt"]), method = "seasonal",
                  T_period = 52.18,
                  num_surr = num_surr,
                  alpha=alpha) %>%
    as.data.frame()
  df=df %>% select(-"plt")
  surrA=bind_cols(df,surr_data)
}


################################### Determine optimal E for the system by each state ###########################################

jishu <- function(x){
  ifelse(x%%2 ==0,F,T)
}

make_pred_nozeroL <- function(dat){
  dat <- dat %>% mutate(year=year(date),month=month(date),day=day(date),n=1:nrow(dat))
  dat1 <- dat %>% filter(month==5) %>% group_by(year) %>% filter(day==max(day))
  dat2 <- dat %>% filter(month==10) %>% group_by(year) %>% filter(day==min(day))
  I_zero_strings <- c(dat1$n,dat2$n)[order(c(dat1$n,dat2$n))]
  
  if(I_zero_strings[1]!=1 | jishu(length(I_zero_strings))==T) {
    I_zero_strings=c(1,I_zero_strings)
  } else {
    I_zero_strings=I_zero_strings
  }
  
  if (I_zero_strings[2]==1) {
    I_zero_strings=I_zero_strings[3:length(I_zero_strings)]
  }
  
  lib_out <- matrix(I_zero_strings,ncol=2,byrow=T)
  return(lib_out)
}

fn_E_smapc=function(data,ST,y){
  
  df=data %>%
    filter(state==ST) %>%
    select(date,y) %>%
    na.omit()
  
  M <- nrow(df)
  
  lib <- make_pred_nozeroL(df)
  pred <- make_pred_nozeroL(df)
  
  E=EmbedDimension(dataFrame=df,
                   lib=lib,
                   pred=pred,
                   columns=y, target=y,
                   maxE=20,
                   showPlot=F)
  temp=data.frame(dis=y,E,ST)
  temp}

plist=list(data=list(df_smapc),y=c('fluP'),ST=unique(dfA$state)) %>%
  cross_df()

E_smapc_out=plist %>% pmap_df(fn_E_smapc)

E_smapc =E_smapc_out %>% filter(E %in% 2:6) %>%
  group_by(dis,ST) %>%
  filter(rho==max(rho))  %>%
  as.data.frame() 

################################### Determine optimal theta for S-map by each state ###########################################

fn_theta_justY=function(data, ST, dis, theta){
  
  E=E_smapc[E_smapc[,'ST']==ST & E_smapc[,'dis']==dis,'E']
  
  df=data %>%
    filter(state==ST) %>%
    select(date,dis)
  
  M <- nrow(df)
  
  lib=make_pred_nozeroL(df)
  pred=make_pred_nozeroL(df)
  
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


############################################# surrugate CCM #####################################################
fn_season_ccm=function(data,ST,x,y,tp_value){
  
  df=data %>%
    filter(state==ST) %>%
    select(date,y,x)
  
  E=E_smapc[E_smapc[,'ST']==ST & E_smapc[,'dis']==y,'E']
  
  alpha=sd_data %>% filter(states==ST & plt==x) %>%
    select(sd_PNAS) %>% pull()
  
  surr_data <- fn_surr_data(dfB,ST,x,0)
  
  all_data <- df %>% left_join(surr_data,by="date")
  
  names(all_data) = c("date", y, 'T1',paste0("T", 2:(num_surr+1)))
  
  m=nrow(all_data) %>% as.data.frame()
  
  libSize =c(E+2,m-E-2)
  
  rho_surr <- NULL
  
  for (i in 1:(num_surr+1)) {
    targetCol = paste("T", i, sep = "")
    ccm_out = CCM(dataFrame = all_data, E = E, Tp = tp_value,
                  columns = y,
                  target = targetCol,
                  libSizes = libSize,
                  random=T,
                  sample = num_sample,
                  seed=2019)
    col = paste(y, ":", targetCol, sep = "")
    dat=ccm_out %>% select(LibSize,col)
    names(dat)=c("lib","rho")
    test1=mutate(dat,i=i,dis=y,plt=x,E=E,
                 tp_value=tp_value,state=ST)
    rho_surr <- rbind(rho_surr,test1)
    
  }
  rho_surr
}

plist=list(data=list(dfA),
           ST=unique(dfA$state),
           y=c('fluP'),
           x=c('o3',"temp","ah"),
           tp_value=-2:0) %>% cross_df()

cl <- makeCluster(cores[1],type = core_type)
registerDoParallel(cl)
ccm_out=foreach(j = 1:nrow(plist),
                .packages = c("rEDM","tidyverse"),
                .combine=rbind,
                # .export='num_sample',
                .inorder=FALSE)  %dopar% {
                  ccm_out=plist[j,] %>% pmap_df(fn_season_ccm)
                }
stopCluster(cl)

dat_min <- ccm_out %>% filter(lib<50)
dat_min <- dat_min[order(dat_min$state,dat_min$tp_value,dat_min$plt,dat_min$dis,dat_min$i),]
dat_max <- ccm_out %>% filter(lib>50)
dat_max <- dat_max[order(dat_max$state,dat_max$tp_value,dat_max$plt,dat_max$dis,dat_max$i),]

dat <- cbind(dat_max,dat_min[,"rho"])

names(dat) <- c("lib","rho_max","i","dis","plt","E","tp_value","ST","rho_min")

dat1 <- dat %>% mutate(rho=rho_max-rho_min)
ccm_out <- dat1

ccm_out_raw = ccm_out  %>%
  filter(i==1) %>%
  select(plt, ST,tp_value, rho)

ccm_p=ccm_out %>%
  group_by(dis,plt,E,ST,tp_value) %>%
  summarise(p=1-ecdf(rho[i != 1])(rho[i == 1])) %>%
  left_join(ccm_out_raw, by=c("plt", "ST", "tp_value")) %>%
  rename(rho_raw=rho) %>%
  arrange(dis,plt)

ccm_p %>% filter(p==0)
ccm_p %>% filter(rho_raw<=0)

ccm_p=ccm_p %>%
  mutate(p=ifelse(p==0, 0.0005, p)) %>%
  mutate(p=ifelse(rho_raw<=0, 1, p))

ccm_p_table=ccm_p %>% 
  select(ST, Pollutant=plt, Lag=tp_value, CCM_skill=rho_raw, P=p) %>% 
  mutate(Pollutant=factor(Pollutant, levels=c("o3", "ah", "temp"), labels=c('O3', 'AH', 'T')))

knitr::kable(ccm_p_table, caption = "CCM causality test results of Kentucky.")

################################################ Effect strength #################################################
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
  
  lib=make_pred_nozeroL(dataFrame)
  pred=make_pred_nozeroL(dataFrame)
  
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

tempA <- C_out %>%   group_by(tp_value,E,dis,plt,ST) %>%
  summarise(effect=mean(effect,na.rm=T))

tempA <- tempA %>% filter(plt=="o3") %>% mutate(ST="KY")

ggplot(tempA,aes(tp_value,effect))+geom_line()+geom_hline(yintercept = 0) +
  ylab("Effect size")+
  xlab("lag")+
  facet_wrap(~ST,scales='free',nrow = 6)+
  theme_bw() +
  theme(axis.text.x=element_text(size = 12, color= "black"),
        axis.text.y=element_text(size = 12, color= "black"),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        axis.ticks.length=unit(0.2,'cm'),
        strip.text.x = element_text(size = 15),
        legend.position="none",
        panel.background = element_blank(),
        plot.title = element_text(size = 18, hjust=0.5))
