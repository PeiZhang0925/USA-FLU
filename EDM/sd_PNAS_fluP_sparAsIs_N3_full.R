rm(list=ls())


#+ packages
packages=c('tidyverse','knitr','lubridate','rEDM',
           'tsModel','doParallel','foreach')
lapply(packages, require, character.only=T)

##########################################################
load('usa_Flu_N3_proxy_Data_full.rda')

dfA=usa_Flu_N3_proxy_Data_full

df_flu=dfA %>% select(date,state,ah,o3,temp) %>%
  group_by(state) %>%
  gather(key,value,-c(date,state))

ggplot(df_flu,aes(date,value))+geom_line()+
  facet_grid(key~state,scales='free')



##### To estimate Spar first ############

# https://stackoverflow.com/questions/14929268/how-do-i-select-the-smoothing-parameter-for-smooth-spline




########################################################

cores_all=detectCores()
cores=ifelse(cores_all<9,4,10)

# core_type='FORK'
core_type='PSOCK'
########################################################


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
    select(value) %>% 
    as.vector() %>% pull() 
  
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
  
  # plot(spars, ss, 'l', xlab = 'spar', ylab = 'Cross Validation Residual Sum of Squares' , main = 'CV RSS vs Spar')
  spar=spars[which.min(ss)]
  
  data.frame(state=states,plt=keys,spar)
}


########################################################
plist=list(states=unique(df_flu$state),
           keys=c('o3','ah','temp')) %>% 
  cross_df()

cl <- makeCluster(cores[1],type = core_type) 
registerDoParallel(cl)

flu_spar=plist %>% pmap_df(fn_state_spar)

stopCluster(cl)
########################################################

flu_spar 

###### PNAS function #########
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
########################################################

fn_anomaly_PNAS=function(states,plts){
  vec_t=df_flu %>% 
    ungroup() %>% 
    filter(state==states & key==plts) %>%
    select(date) %>% as.vector() %>% pull()
  
  vec_x=df_flu %>% 
    ungroup() %>% 
    filter(state==states & key==plts) %>%
    select(value) %>% 
    as.vector() %>% pull() 
  
  spars=flu_spar %>% filter(state==states & plt==plts) %>%
    select(spar) %>%
    as.vector() %>% pull()
  # spars=0.8
  
  df_9=yearday_anom(vec_t,vec_x,spars) 
  sd=sd(df_9$anomaly,na.rm=TRUE)
  data.frame(states,plt=plts,spar=spars,sd_PNAS=sd)
}


plist=list(states=unique(df_flu$state),
           plts=c('o3','ah','temp')) %>% 
  cross_df()


sd_PNAS=plist %>% pmap_df(fn_anomaly_PNAS)

sd_PNAS_fluP_sparAsIs_Not08=sd_PNAS 

sd_PNAS_fluP_sparAsIs_N3_full= sd_PNAS_fluP_sparAsIs_Not08

save(sd_PNAS_fluP_sparAsIs_N3_full,file='sd_PNAS_fluP_sparAsIs_N3_full.rda')

#' sd_mean calculation was given up because there are not 
#' sufficient number of identical 'days' for averaging 
#' 


