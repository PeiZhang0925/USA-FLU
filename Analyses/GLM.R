# ----------------------
# 1. Load required packages
# ----------------------
library(tidyverse)
library(lubridate)
library(mgcv)
library(tsModel)

# ----------------------
# 2. Load data
# ----------------------
df <- read.csv("state_ky.csv")

dfA <- df %>% 
  mutate(
    time = as.numeric(date),
    year = year(date),
    month = month(date)
  )

# ----------------------
# 3. Define GAM functions
# ----------------------
gam_function <- function(dat, Y='fluP', x1='value', xlag=0, cityname, variable) {
  
  df <- dat %>%
    filter(state == cityname) %>%
    mutate(
      loglag1 = log(Lag(fluP, 1)),
      loglag2 = Lag(loglag1, 1)
    )
  
  formula_str <- switch(variable,
                        o3   = paste(Y, '~Lag(', x1, ',', xlag, ')+ Lag(ah,', xlag, ')+ Lag(temp,', xlag, ')+ as.factor(year) + as.factor(month) + loglag1', sep=''),
                        ah   = paste(Y, '~Lag(', x1, ',', xlag, ')+ Lag(o3,', xlag, ')+ Lag(temp,', xlag, ')+ as.factor(year) + as.factor(month) + loglag1', sep=''),
                        temp = paste(Y, '~Lag(', x1, ',', xlag, ')+ Lag(o3,', xlag, ')+ Lag(ah,', xlag, ')+ as.factor(year) + as.factor(month) + loglag1', sep='')
  )
  
  fit <- gam(as.formula(formula_str), data=df, family=quasibinomial, na.action=na.omit)
  summ <- summary(fit)
  
  outA <- data.frame(
    beta = summ$p.coeff[2],
    se   = summ$se[2],
    t    = summ$p.t[2],
    p    = summ$p.pv[2],
    state = cityname,
    lag   = xlag
  )
  
  return(outA)
}

# ----------------------
# 4. GAM analysis
# ----------------------
variables <- c("o3", "ah", "temp")
results <- list()

for (var in variables) {
  plist <- list(
    dat      = list(dfA),
    Y        = 'fluP',
    x1       = var,
    xlag     = 0:2,
    cityname = unique(dfA$state),
    variable = var
  ) %>% cross_df()
  
  results[[var]] <- plist %>% pmap_df(gam_function)
}

# Store the results for easy access
out_o3   <- results$o3
out_ah   <- results$ah
out_temp <- results$temp

# calculate overall and state-specific SD 

df_SD_ST <- dfA %>% group_by(state) %>%
  dplyr::summarize(o3=sd(o3),
                   ah=sd(ah),
                   temp=sd(temp)) %>%
  gather(plt, SD, -state) %>%
  mutate(state=toupper(state))

# transform beta to effect size of each SD change
df1=out_o3  %>% as.data.frame() %>% mutate(plt='o3')
df2=out_ah  %>% as.data.frame() %>% mutate(plt='ah')
df3=out_temp %>% as.data.frame() %>% mutate(plt='temp')

df_gam_final=rbind(df1,df2,df3) %>%
  mutate(state = toupper(state)) %>%
  left_join(df_SD_ST,by=c("state",'plt')) %>%
  mutate(betalow=beta-1.96*se, betahigh=beta+1.96*se) %>%
  mutate(Size=beta*SD,SizeL=betalow*SD,SizeH=betahigh*SD) %>%
  select(state, plt, lag, 
         Size, SizeL, SizeH,  gam.p=p) %>%
  mutate(sig=ifelse(gam.p<0.05, 'sig', 'non_sig'),
         sig=factor(sig, levels=c("non_sig","sig")))

