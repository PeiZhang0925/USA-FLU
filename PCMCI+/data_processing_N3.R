rm(list = ls())

#+ packages, include=FALSE 
packages=c('tidyverse','knitr','lubridate')
lapply(packages, require, character.only=T) 
#- 


load('usa_Flu_N3_proxy_Data.rda')
load('usa_Flu_N3_proxy_Data_full.rda')

df1 <- usa_Flu_N3_proxy_Data_full %>% filter(state %in% c("ak","hi"))
df2 <- usa_Flu_N3_proxy_Data %>% filter(!state %in% c("ak","hi"))


df <- rbind(df1,df2[,1:12])

logitTransform <- function(p) { log(p/(1-p)) }

dfA=df %>%
  mutate(logitfluP=logitTransform(fluP+0.0000001)) 
dfA$fluP <- dfA$logitfluP

dfA=dfA %>% 
  mutate(month=month(date)) %>% 
  # mutate(season=as.numeric(month %in% c(3:11))) %>% 
  mutate(season=
           ifelse((state=='hi' & month %in% c(5:10)) | 
                    (state=='ak' & month %in% c(5:9)) | 
                    (state !='hi' & state !='ak' & month %in% c(3:11)), 
                  1,0) )

table(dfA$state)
# dfA <- dfA %>% filter(!month %in% c(6:9))
###################

fn_norm=function(x){
  x=as.numeric(x)
  xt=x
  # Normalization (zero mean & unity variance)
  xt=(xt-mean(xt,na.rm=T))/sd(xt,na.rm=T)
  return(xt)}
###################
###################

new.row <- data.frame(
  date=rep(as.Date(c('2015-10-04')),45),
  state=rep(unique(dfA$state),each=1),season=999.,
  ah=999.,o3=999.,temp=999.,fluP=999.)

new.row1 <- data.frame(
  date=rep(as.Date(c('2011-07-04')),45),
  state=rep(unique(dfA$state),each=1),season=999.,
  ah=999.,o3=999.,temp=999.,fluP=999.)

new.row2 <- data.frame(
  date=rep(as.Date(c('2012-07-04')),45),
  state=rep(unique(dfA$state),each=1),season=999.,
  ah=999.,o3=999.,temp=999.,fluP=999.)

new.row3 <- data.frame(
  date=rep(as.Date(c('2013-07-04')),45),
  state=rep(unique(dfA$state),each=1),season=999.,
  ah=999.,o3=999.,temp=999.,fluP=999.)

new.row4 <- data.frame(
  date=rep(as.Date(c('2014-07-04')),45),
  state=rep(unique(dfA$state),each=1),season=999.,
  ah=999.,o3=999.,temp=999.,fluP=999.)

new.row5 <- data.frame(
  date=rep(as.Date(c('2015-07-04')),45),
  state=rep(unique(dfA$state),each=1),season=999.,
  ah=999.,o3=999.,temp=999.,fluP=999.)

###########################


usaflu_raw=dfA %>% 
  select(date,state,season, ah,o3,temp,fluP) %>% 
  bind_rows(new.row1) %>%
  bind_rows(new.row2) %>%
  bind_rows(new.row3) %>%
  bind_rows(new.row4) %>%
  bind_rows(new.row5) %>%
  arrange(state,date) %>% 
  rownames_to_column("row")
table(usaflu_raw$state)

usaflu_norm=dfA %>% 
  select(date,state,season, ah,o3,temp,fluP) %>%
  group_by(state) %>% 
  mutate_at(vars(-date,-state,-season),fn_norm) %>%
  bind_rows(new.row) %>%
  bind_rows(new.row1) %>%  
  bind_rows(new.row2) %>%
  bind_rows(new.row3) %>%
  bind_rows(new.row4) %>%
  bind_rows(new.row5) %>%
  arrange(state,date) %>% 
  arrange(state,date) %>% 
  rownames_to_column("row") %>% 
  filter(!state %in% c("hi","ak" ) ) 


usaflu_raw <- usaflu_raw %>% mutate(state=case_when(state == "ak" ~ "Alaska", state == "az" ~ "Arizona",
                                                    state == "ca" ~ "California", state == "co" ~ "Colorado", 
                                                    state == "ga" ~ "Georgia", state == "hi" ~ "Hawaii",
                                                    state == "sd" ~ "South Dakota", state == "il" ~ "Illinois",
                                                    state == "ky" ~ "Kentucky", state == "mn" ~ "Minnesota", 
                                                    state == "mo" ~ "Missouri", state == "ny" ~ "New York",
                                                    state == "oh" ~ "Ohio", state == "pa" ~ "Pennsylvania",
                                                    state == "tx" ~ "Texas", state == "ut" ~ "Utah", 
                                                    state == "va" ~ "Virginia", state == "wa" ~ "Washington",
                                                    state == "wi" ~ "Wisconsin", state == "wv" ~ "West Virginia",
                                                    state == "ar" ~ "Arkansas", state == "la" ~ "Louisiana", 
                                                    state == "sc" ~ "South Carolina", state == "de" ~ "Delaware",
                                                    state == "ct" ~ "Connecticut", state == "sd" ~ "South Dakota",
                                                    state == "al" ~ "Alabama", state == "md" ~ "Maryland",
                                                    state == "ma" ~ "Massachusetts", state == "or" ~ "Oregon",
                                                    state == "ia" ~ "Iowa", state == "mt" ~ "Montana",
                                                    state == "nm" ~ "New Mexico", state == "tn" ~ "Tennessee",
                                                    state == "ok" ~ "Oklahoma", state == "nv" ~ "Nevada",
                                                    state == "nc" ~ "North Carolina", state == "ms" ~ "Mississippi",
                                                    state == "wy" ~ "Wyoming", state == "me" ~ "Maine",
                                                    state == "nd" ~ "North Dakota", state == "nh" ~ "New Hampshire",
                                                    state == "id" ~ "Idaho", state == "ks" ~ "Kansas",
                                                    state == "vt" ~ "Vermont", state == "mi" ~ "Michigan"
))


write_csv(usaflu_raw, "usaflu_raw_N3_logit_fluseason.csv")
write_csv(usaflu_norm, "usaflu_norm_N3_logit_fluseason.csv")
