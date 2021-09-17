rm(list = ls())

packages=c('tidyverse','knitr','lubridate')
lapply(packages, require, character.only=T) 

load('usa_Flu_J_proxy_Data.rda')

dfA=usa_Flu_J_proxy_Data %>% 
  mutate(month=month(date)) %>% 
  mutate(season=
           ifelse(state=='hi' & month %in% c(5:10),1,
                  as.numeric(month %in% c(6:8))))

###################
fn_norm=function(x){
  x=as.numeric(x)
  xt=x
  # Normalization (zero mean & unity variance)
  xt=(xt-mean(xt,na.rm=T))/sd(xt,na.rm=T)
  return(xt)}

new.row <- data.frame(
  date=rep(as.Date(c('2015-10-04')),30),
  state=rep(unique(dfA$state),each=1),
  ah=999.,o3=999.,temp=999.,fluP=999.)

###########################
usaflu_raw=dfA %>% 
  select(date,state,season, ah,o3,temp,fluP) %>%
  group_by(state) %>% 
  arrange(state,date) %>% 
  rownames_to_column("row")


usaflu_norm=dfA %>% 
  select(date,state,season, ah,o3,temp,fluP) %>%
  group_by(state) %>% 
  mutate_at(vars(-date,-state),fn_norm) %>%
  bind_rows(new.row) %>%
  arrange(state,date) %>% 
  rownames_to_column("row")


write_csv(usaflu_raw, "usaflu_raw.csv")
write_csv(usaflu_norm, "usaflu_norm.csv")
