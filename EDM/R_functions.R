logitTransform <- function(p) { log(p/(1-p)) }

normFunc=function(x){(x-mean(x, na.rm = T))/sd(x, na.rm = T)}

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
  
  # I_zero_strings=I_zero_strings[1:(length(I_zero_strings)-1)]
  lib_out <- matrix(I_zero_strings,ncol=2,byrow=T)
  return(lib_out)
}

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
  
  # plot(spars, ss, 'l', xlab = 'spar', ylab = 'Cross Validation Residual Sum of Squares' , main = 'CV RSS vs Spar')
  spar=spars[which.min(ss)]
  
  data.frame(state=states,plt=keys,spar)
}

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
  # spars=0.8
  
  df_9=yearday_anom(vec_t,vec_x,spars)
  sd=sd(df_9$anomaly,na.rm=TRUE)
  data.frame(states,plt=plts,spar=spars,sd_PNAS=sd)
}

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
                      embed_1[E:M, 1:(E-1)],df[E:M, 'plt_1']  ) %>%
    as.data.frame()
  
  names(dataFrame_1)=c('date',dis,letters[1:(E-1)],plt_1)
  columns_1 = paste(paste(letters[1:(E-1)],collapse =' '), plt_1 ,
                    sep=' ')
  
  m = nrow(dataFrame_1)
  
  theta=best_theta[best_theta$state==ST & best_theta$dis==dis,"theta"]
  
  smap_1 = SMap(dataFrame = dataFrame_1,
                embedded = TRUE,
                columns = columns_1,
                target = dis,
                lib = make_pred_nozeroL(dataFrame_1),
                pred = make_pred_nozeroL(dataFrame_1),
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
                      embed_1[E:M, 1:(E-1)],df[E:M, 'plt_1'],
                      df[E:M, 'plt_2']) %>% as.data.frame()
  
  names(dataFrame_2)=c('date',dis,letters[1:(E-1)],plt_1,plt_2)
  columns_2 = paste(paste(letters[1:(E-1)],collapse =' '), plt_1,plt_2 ,
                    sep=' ')
  
  m = nrow(dataFrame_2)
  
  theta=best_theta[best_theta$state==ST & best_theta$dis==dis,"theta"]
  
  smap_2 = SMap(dataFrame = dataFrame_2,
                embedded = TRUE,
                columns = columns_2,
                target = dis,
                lib = make_pred_nozeroL(dataFrame_2),
                pred = make_pred_nozeroL(dataFrame_2),
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

one_sample_mean=function(x){
  #x=df_123_fluP %>% filter(tp_values==tp_value, plt == plts) %>% select(delta_rho) %>% pull()
  wilcox.test(x, mu = 0, alternative = "greater")
}

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
    
    surr_data <- fn_surr_data(dfB_smapc,ST,plt_1,tp_value)
    
    all_data <- df %>% left_join(surr_data,by="date")
    
    names(all_data) <- c("date", dis, 'T1',paste0("T", 2:(num_surr+1)))
    M <- nrow(all_data)
    out <- foreach(i = 1:(num_surr+1),
                   .packages = c("rEDM","tidyverse",'lubridate'),
                   .export=c('num_sample','num_surr','best_theta'),
                   .combine=rbind,
                   .inorder=FALSE)  %dopar% {
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
                     targetCol <- paste("T", i, sep = "")
                     embed_1 <- Embed(dataFrame = all_data, E = E, tau = -1, columns = dis )
                     dataFrame_1 <- cbind(all_data[E:M, 'date'],all_data[E:M, dis],embed_1[E:M, 1:(E-1)],
                                          all_data[E:M, targetCol]) %>% as.data.frame()
                     
                     # dataFrame_1 <- cbind(all_data[E:M, 'date'],all_data[E:M, dis],embed_1[E:M, 2:(E-0)],
                     #                      all_data[E:M, targetCol]) %>% as.data.frame()
                     
                     names(dataFrame_1) <- c('date',dis,letters[1:(E-1)],plt_1)
                     columns_1 <- paste(paste(letters[1:(E-1)],collapse =' '), plt_1 , sep=' ')
                     theta=best_theta[best_theta$state==ST & best_theta$dis==dis,"theta"]
                     smap_1 <- SMap(dataFrame = dataFrame_1,embedded = TRUE,
                                    columns = columns_1,target = dis,
                                    lib = make_pred_nozeroL(dataFrame_1),
                                    # pred = c(1, nrow(dataFrame_1)),
                                    pred = make_pred_nozeroL(dataFrame_1),
                                    theta = theta,Tp = 1)
                     out_1=compute_stats(smap_1$predictions$Observations,
                                         smap_1$predictions$Predictions) %>% mutate(mfi_012=1)
                     out=rbind(out_1) %>% mutate(i=i,ST=ST,plt=plt_1,dis=dis,tp_value=tp_value,
                                                 alpha=alpha)
                     out
                   }
    out
  }

fn_metap=function(var,tp_values,plts,diss){
  df=filter(var,tp_value==tp_values & plt==plts & dis==diss)
  out=allmetap(df$p, method = "sumlog") %>% as.data.frame()
  mutate(out,tp_value=tp_values,plt=plts,dis=diss)
}

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

p_function=function(x,pH,pL){
  dfg2_Q1=dfg2 %>% filter(grp=='surr') %>%
    group_by(dis,plt,tp_value,E) %>%
    summarise(ph=quantile(.data[[x]], pH,na.rm=T),
              pl=quantile(.data[[x]], pL,na.rm=T))
  
  dfg2_Q0=dfg2 %>% filter(grp=='raw') %>% select(-grp) %>%
    group_by(dis,plt,tp_value,E) %>%
    select(p0=x)
  
  dfg2_Q=left_join(dfg2_Q1,dfg2_Q0,by=c('E','tp_value',
                                        'dis','plt'))
  
  ggplot(dfg2_Q) +
    # geom_hline(yintercept = 0,linetype='dashed') +
    geom_line(aes(tp_value,p0),size=0.5)+
    geom_point(aes(tp_value,p0),size=1)+
    geom_linerange(aes(x=tp_value,ymin =pl, ymax = ph),
                   size=1,color='green',alpha=0.5)+
    theme_bw() +
    ylab('Forecast improvement')+
    xlab("Lag")+
    theme(axis.text.x=element_text(size = 12, color= "black"),
          axis.text.y=element_text(size = 12, color= "black"),
          axis.title.x=element_text(size=15),
          axis.title.y=element_text(size=15),
          axis.ticks.length=unit(0.2,'cm'),
          legend.position="none",
          panel.background = element_blank(),
          plot.title = element_text(size = 18, hjust=0.5)) +
    facet_grid(plt~dis,scales='free')
}

