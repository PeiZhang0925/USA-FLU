
rm(list = ls())

#+ packages, include=FALSE 
packages=c('tidyverse','knitr','lubridate')
lapply(packages, require, character.only=T) 

load("Data_fluseason.rda")

# logit transformation
dfA=usa_Flu_P_proxy_Data_B 

logitTransform <- function(p) { log(p/(1-p)) }

dfA$logitfluP <- logitTransform(dfA$fluP)
dfA$fluP <- dfA$logitfluP

# centralization
siteminus=function(x){x=x-mean(x,na.rm=T)}

statemean=function(x){x=mean(x,na.rm=T)}

dfA=dfA %>% 
  select(date,state,ah,o3,temp,fluP) %>% 
  pivot_longer(cols=(ah:temp)) %>% 
  group_by(state,name) %>%
  mutate(siteminus=siteminus(value)) %>%
  group_by(name) %>%
  mutate(statemean=statemean(value)) %>% 
  mutate(value_centr=siteminus+statemean) %>% 
  select(-c(value,siteminus,statemean)) %>% 
  pivot_wider(names_from='name',values_from='value_centr') 

# add season variable
dfA=dfA %>%
  mutate(month=month(date)) %>%
  mutate(season=
           ifelse((state=='hi' & month %in% c(5:10)) |
                    (state=='ak' & month %in% c(5:9)) |
                    (state !='hi' & state !='ak' & month %in% c(3:11)),
                  1,0) )

# add NA lines to separate discontinuous data (different states or discontinuous in date)
new.row1 <- data.frame(
  date=rep(as.Date(c('2011-07-04')),46),
  state=rep(unique(dfA$state),each=1),season=999.,
  ah=999.,o3=999.,temp=999.,fluP=999.)

new.row2 <- data.frame(
  date=rep(as.Date(c('2012-07-04')),46),
  state=rep(unique(dfA$state),each=1),season=999.,
  ah=999.,o3=999.,temp=999.,fluP=999.)

new.row3 <- data.frame(
  date=rep(as.Date(c('2013-07-04')),46),
  state=rep(unique(dfA$state),each=1),season=999.,
  ah=999.,o3=999.,temp=999.,fluP=999.)

new.row4 <- data.frame(
  date=rep(as.Date(c('2014-07-04')),46),
  state=rep(unique(dfA$state),each=1),season=999.,
  ah=999.,o3=999.,temp=999.,fluP=999.)

new.row5 <- data.frame(
  date=rep(as.Date(c('2010-10-01')),46),
  state=rep(unique(dfA$state),each=1),season=999.,
  ah=999.,o3=999.,temp=999.,fluP=999.)

# state-wise analysis data
usaflu_raw=dfA %>% 
  select(date,state,season, ah,o3,temp,fluP) %>% 
  bind_rows(new.row1) %>%
  bind_rows(new.row2) %>%
  bind_rows(new.row3) %>%
  bind_rows(new.row4) %>%
  bind_rows(new.row5) %>%
  arrange(state,date) %>% 
  rownames_to_column("row")

# overall analysis data (combined all the states and normalized)
fn_norm=function(x){
  x=as.numeric(x)
  xt=x
  # Normalization (zero mean & unity variance)
  xt=(xt-mean(xt,na.rm=T))/sd(xt,na.rm=T)
  return(xt)}

usaflu_norm=dfA %>% 
  select(date,state,season, ah,o3,temp,fluP) %>%
  ungroup() %>% 
  mutate_at(vars(-date,-state,-season),fn_norm) %>%
  bind_rows(new.row1) %>%  
  bind_rows(new.row2) %>%
  bind_rows(new.row3) %>%
  bind_rows(new.row4) %>%
  bind_rows(new.row5) %>%
  arrange(state,date) %>% 
  rownames_to_column("row")


usaflu_raw <- usaflu_raw %>% mutate(state=case_when(state == "ak" ~ "AK", state == "al" ~ "AL",
                                                    state == "ar" ~ "AR",state == "az" ~ "AZ",
                                                    state == "ca" ~ "CA", state == "co" ~ "CO",
                                                    state == "ct" ~ "CT", state == "de" ~ "DE",
                                                    state == "ga" ~ "GA", state == "hi" ~ "HI",
                                                    state == "ia" ~ "IA",state == "id" ~ "ID",
                                                    state == "il" ~ "IL",state == "in" ~ "IN",
                                                    state == "ks" ~ "KS",state == "ky" ~ "KY",
                                                    state == "la" ~ "LA", state == "ma" ~ "MA", 
                                                    state == "md" ~ "MD",state == "me" ~ "ME",
                                                    state == "mi" ~ "MI",state == "mn" ~ "MN", 
                                                    state == "mo" ~ "MO",state == "ms" ~ "MS",
                                                    state == "mt" ~ "MT",state == "nc" ~ "NC", 
                                                    state == "nd" ~ "ND",state == "ne" ~ "NE",
                                                    state == "nh" ~ "NH",state == "nm" ~ "NM",
                                                    state == "nv" ~ "NV", state == "ny" ~ "NY",
                                                    state == "oh" ~ "OH",state == "ok" ~ "OK",
                                                    state == "or" ~ "OR",state == "pa" ~ "PA",
                                                    state == "sc" ~ "SC", state == "sd" ~ "SD",
                                                    state == "tn" ~ "TN",state == "tx" ~ "TX",
                                                    state == "ut" ~ "UT",state == "va" ~ "VA",
                                                    state == "wa" ~ "WA",state == "wi" ~ "WI",
                                                    state == "wv" ~ "WV", state == "wy" ~ "WY"
                                                    
))

write_csv(usaflu_raw, "usaflu_raw.csv")
write_csv(usaflu_norm, "usaflu_norm.csv")
