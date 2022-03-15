################################## Load packages ###################################
options(warn = -1)
options(width=80) 
options(stringsAsFactors=FALSE)
options(scipen = 6) # bias against scientific notation
options(digits = 3) # show fewer decimal places

library(tidyverse); library(lubridate); library(dlnm); library(splines);
library(tsModel); library(gnm);library(ggpubr);library(metafor);
library(mgcv)

################################## Load data ###################################
load("usa_Flu_N3_proxy_Data.rda")

dfA=usa_Flu_N3_proxy_Data
dfA <- dfA %>% 
  mutate(time=as.numeric(date),
         year=year(date), 
         month=month(date),
         week=week(date),
         ym=paste(year(date), month(date),sep = "")) %>%
  mutate(season=ifelse(state=='hi' & month %in% c(5:10) | 
                         state=='ak' & month %in% c(5:9) | 
                         state !='hi' & state !='ak' & month %in% c(3:11), 
                       1,0 )) %>%
  rename(ST=state) %>% filter(!ST %in% c("hi","ak" ) )

dfA$state <- dfA$ST
  
# dfA <- dfA %>% mutate(state=case_when(state == "ak" ~ "Alaska", state == "az" ~ "Arizona",
#                                                       state == "ca" ~ "California", state == "co" ~ "Colorado", 
#                                                       state == "ga" ~ "Georgia", state == "hi" ~ "Hawaii",
#                                                       state == "sd" ~ "South Dakota", state == "il" ~ "Illinois",
#                                                       state == "ky" ~ "Kentucky", state == "mn" ~ "Minnesota", 
#                                                       state == "mo" ~ "Missouri", state == "ny" ~ "New York",
#                                                       state == "oh" ~ "Ohio", state == "pa" ~ "Pennsylvania",
#                                                       state == "tx" ~ "Texas", state == "ut" ~ "Utah", 
#                                                       state == "va" ~ "Virginia", state == "wa" ~ "Washington",
#                                                       state == "wi" ~ "Wisconsin", state == "wv" ~ "West Virginia",
#                                                       state == "ar" ~ "Arkansas", state == "la" ~ "Louisiana", 
#                                                       state == "sc" ~ "South Carolina", state == "de" ~ "Delaware",
#                                                       state == "ct" ~ "Connecticut", state == "sd" ~ "South Dakota",
#                                                       state == "al" ~ "Alabama", state == "md" ~ "Maryland",
#                                                       state == "ma" ~ "Massachusetts", state == "or" ~ "Oregon",
#                                                       state == "ia" ~ "Iowa", state == "mt" ~ "Montana",
#                                                       state == "nm" ~ "New Mexico", state == "tn" ~ "Tennessee",
#                                                       state == "ok" ~ "Oklahoma", state == "nv" ~ "Nevada",
#                                                       state == "nc" ~ "North Carolina", state == "ms" ~ "Mississippi",
#                                                       state == "wy" ~ "Wyoming", state == "me" ~ "Maine",
#                                                       state == "nd" ~ "North Dakota", state == "nh" ~ "New Hampshire",
#                                                       state == "id" ~ "Idaho", state == "ks" ~ "Kansas",
#                                                       state == "vt" ~ "Vermont", state == "mi" ~ "Michigan"
#   ))

dfA$ST <- dfA$state
dfA <- dfA %>% filter(!month %in% c(6:9))
################################## State-specific GAM analysis ###################################

fgam2=function(Y='cvsd',x1='value',x2,x3,xlag=0,cityname){
  
  df <- dfA %>% filter(state == cityname) %>% 
    mutate(loglag1 = log(Lag(fluP, 1)+0.0000001)) 
  
  # s(time,k=(5*8+1),fx=T), as.factor(year) + as.factor(month)
  
  fit=gam(as.formula(paste(Y,
                           '~Lag(',x1,',',xlag,')+as.factor(year) + as.factor(month)+loglag1', sep='')),
          data=df,family=quasibinomial,na.action=na.omit, 
          #weights=w,
          control = glm.control(epsilon = 1e-10, maxit = 1000))
  
  summ=summary(fit)
  
  outA=data.frame(summ$p.coeff[2],summ$se[2],summ$p.t[2],summ$p.pv[2])
  names(outA)=c('beta','se','t','p')
  outA$state=cityname
  outA$lag=xlag
  outA$plt=x1
  outA
}

plist=list(Y='fluP',x1=c('o3','temp','ah'),x2='temp',x3='ah',
           xlag=0:2,
           cityname=unique(dfA$state)) %>% cross_df()
out <- plist %>% pmap_df(fgam2)

out

GAMres <- out %>% select(state,lag, effect=beta, p,plt) %>% as.tibble()


################################## Meta analysis of state-specific results ###################################

result <- matrix(nrow=3, ncol=7, dimnames=list(1:3, c("state", "lag","effect",'p','lb','ub',"plt")))

for (j in c('o3')) {
  for (i in 1:3) {
    out.lag <- out %>% filter(!state %in% c("Hawaii","Alaska"))%>% filter(lag==i-1) %>% filter(plt==j)
    meta <- rma(yi=beta, sei=se, slab=state, method="REML", data=out.lag)
    meta.re <- with(meta, c(b, pval, ci.lb,ci.ub))
    result[i,1] = "OVERALL"
    result[i,2] = i-1
    result[i,3] = meta.re[1]
    result[i,4] = meta.re[2]
    result[i,5] = meta.re[3]
    result[i,6] = meta.re[4]
    result[i,7] = j
  }
}

result2 <- matrix(nrow=3, ncol=7, dimnames=list(1:3, c("state", "lag","effect",'p','lb','ub',"plt")))
for (j in c('temp')) {
  for (i in 1:3) {
    out.lag <- out %>% filter(!state %in% c("Hawaii","Alaska"))%>% filter(lag==i-1) %>% filter(plt==j)
    meta <- rma(yi=beta, sei=se, slab=state, method="REML", data=out.lag)
    meta.re <- with(meta, c(b, pval, ci.lb,ci.ub))
    result2[i,1] = "OVERALL"
    result2[i,2] = i-1
    result2[i,3] = meta.re[1]
    result2[i,4] = meta.re[2]
    result2[i,5] = meta.re[3]
    result2[i,6] = meta.re[4]
    result2[i,7] = j
  }
}

result3 <- matrix(nrow=3, ncol=7, dimnames=list(1:3, c("state", "lag","effect",'p','lb','ub',"plt")))
for (j in c('ah')) {
  for (i in 1:3) {
    out.lag <- out %>% filter(!state %in% c("Hawaii","Alaska"))%>% filter(lag==i-1) %>% filter(plt==j)
    meta <- rma(yi=beta, sei=se, slab=state, method="REML", data=out.lag)
    meta.re <- with(meta, c(b, pval, ci.lb,ci.ub))
    result3[i,1] = "OVERALL"
    result3[i,2] = i-1
    result3[i,3] = meta.re[1]
    result3[i,4] = meta.re[2]
    result3[i,5] = meta.re[3]
    result3[i,6] = meta.re[4]
    result3[i,7] = j
  }
}

result <- result %>% as.data.frame() %>% mutate_at(3:6, as.numeric)
result2 <- result2 %>% as.data.frame() %>% mutate_at(3:6, as.numeric)
result3 <- result3 %>% as.data.frame() %>% mutate_at(3:6, as.numeric)
meta.result <- rbind(result,result2,result3)

meta.result

