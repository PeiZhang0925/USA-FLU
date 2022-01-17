rm(list = ls())

#+ packages, include=FALSE 
packages=c('tidyverse','knitr','lubridate')
lapply(packages, require, character.only=T) 
#- 


load('usa_Flu_M_proxy_Data.rda')

dfA=usa_Flu_M_proxy_Data %>% 
  mutate(month=month(date)) %>% 
  # mutate(season=as.numeric(month %in% c(3:11))) %>% 
  mutate(season=
           ifelse((state=='hi' & month %in% c(5:10)) | 
                    (state=='ak' & month %in% c(5:9)) | 
                    (state !='hi' & state !='ak' & month %in% c(3:11)), 
                  1,0) )

table(dfA$state,dfA$season)
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
  date=rep(as.Date(c('2015-10-04')),26),
  state=rep(unique(dfA$state),each=1),season=999.,
  ah=999.,o3=999.,temp=999.,fluP=999.)

###########################


usaflu_raw=dfA %>% 
  select(date,state,season, ah,o3,temp,fluP) %>% 
  filter(!state %in% c("hi","ak" ) ) %>%  
  rownames_to_column("row")
table(usaflu_raw$state)

usaflu_norm=dfA %>% 
  select(date,state,season, ah,o3,temp,fluP) %>%
  group_by(state) %>% 
  mutate_at(vars(-date,-state,-season),fn_norm) %>%
  bind_rows(new.row) %>%
  arrange(state,date) %>% 
  rownames_to_column("row")

usaflu_norm <- usaflu_norm %>% mutate(state=case_when(state == "al" ~ "Alabama", state == "ar" ~ "Arkansas",
                                        state == "az" ~ "Arizona", state == "ca" ~ "California",
                                        state == "co" ~ "Colorado", state == "ct" ~ "Connecticut",
                                         state == "ga" ~ "Georgia", state == "sd" ~ "South Dakota",
                                        state == "il" ~ "Illinois", state == "ky" ~ "Kentucky",
                                        state == "la" ~ "Louisiana", state == "ny" ~ "New York",
                                        state == "mn" ~ "Minnesota", state == "mo" ~ "Missouri",
                                        state == "oh" ~ "Ohio", state == "ok" ~ "Oklahoma",
                                        state == "or" ~ "Oregon", state == "pa" ~ "Pennsylvania",
                                        state == "tx" ~ "Texas", state == "ut" ~ "Utah", 
                                        state == "va" ~ "Virginia", state == "wa" ~ "Washington",
                                        state == "wi" ~ "Wisconsin", state == "wv" ~ "West Virginia",
                                        state == "hi" ~ "Hawaii", state == "ak" ~ "Alaska"))


write_csv(usaflu_raw, "usaflu_raw_M.csv")
write_csv(usaflu_norm, "usaflu_norm_M.csv")
