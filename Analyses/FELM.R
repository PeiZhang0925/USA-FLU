################################## Load packages & data ###################################
library(tidyverse); library(lubridate); library(dlnm); library(splines);
library(tsModel); library(gnm);library(ggpubr);library(metafor);
library(mgcv); library(imputeTS)

library(lfe)

load("Data_fluseason.rda")
df <- usa_Flu_P_proxy_Data_B

holiday <- read.csv("day_public_and_school_holidays_2010_2019.csv")

USholiday = holiday %>% 
  filter(ISO3=="USA") %>% 
  mutate(start_d=cut.Date(as.Date(Date), breaks="weeks", start.on.monday=F)) %>%
  group_by(start_d) %>%
  summarise(all_break=sum(all_break)) %>% 
  mutate(date=as.Date(start_d)) %>% 
  filter(date >= "2010-10-01" & date <= "2015-09-30") %>%
  mutate(holid_ind=ifelse(all_break>2, 1, 0)) # >2 weekend days

dfA <- df %>% 
  mutate(time = as.numeric(date), 
         year = year(date),
         month = month(date),
         season = ifelse(month %in% c(3:5), "spring", ifelse(month %in% c(10:11), "autumn", "winter")),
         loglag1 = log(Lag(fluP, 1)),
         loglag2 = Lag(loglag1, 1)) %>% 
  left_join(USholiday %>% select(date,holid_ind), by="date")

logitTransform <- function(p) { log(p/(1-p)) }
dfA$logitfluP <- logitTransform(dfA$fluP)


################################## State-level FELM ###################################

felm_function_state <- function(dat, Y='logitfluP', cityname, xlag=0) {
  
  df = dat %>% 
    filter(state == cityname) 
  
  formula_str <- paste(Y, '~Lag(o3,', xlag, ')+ Lag(ah,', xlag, ')+ Lag(temp,', xlag, 
                       ')+ as.factor(season) + as.factor(holid_ind) + loglag1 + loglag2 | year', sep='')
  
  fit <- felm(as.formula(formula_str), data=df)
  summ <- summary(fit)
  
  outA <- data.frame(
    plt=c('o3', 'ah', 'temp'),
    beta = summ$coefficients[1:3,1],
    se   = summ$coefficients[1:3,2],
    t    = summ$coefficients[1:3,3],
    p    = summ$coefficients[1:3,4],
    state = cityname,
    lag   = xlag
  )
  
  return(outA)
}

plist <- list(
  dat      = list(dfA),
  Y        = 'logitfluP',
  cityname = unique(dfA$state),
  xlag     = 0:2
) %>% cross_df()

felm_results_state <- plist %>% pmap_df(felm_function_state) %>% 
  as.tibble() %>% 
  mutate(betalow=beta-1.96*se, 
         betahigh=beta+1.96*se) %>% # 95% CI
  mutate(state=as.factor(state),
         lag=factor(lag, levels = c(0,1,2),
                    labels = c("Lag 0","Lag 1","Lag 2")),
         sig=ifelse(p<0.05, 'sig', 'non_sig'),
         sig=factor(sig, levels=c("non_sig","sig"))) %>%
  select(state, plt, lag, beta, betalow, betahigh, p, sig)

################################## nation-level FELM ###################################

felm_function_all <- function(dat, Y='logitfluP', xlag=0) {
    
  formula_str <- paste(Y, '~Lag(o3,', xlag, ')+ Lag(ah,', xlag, ')+ Lag(temp,', xlag, 
                       ')+ as.factor(season) + as.factor(holid_ind) + loglag1 + loglag2 | state + year', sep='')
  
  fit <- felm(as.formula(formula_str), data=dat)
  summ <- summary(fit)
  
  outA <- data.frame(
    plt=c('o3', 'ah', 'temp'),
    beta = summ$coefficients[1:3,1],
    se   = summ$coefficients[1:3,2],
    t    = summ$coefficients[1:3,3],
    p    = summ$coefficients[1:3,4],
    state = "All",
    lag   = xlag
  )
  
  return(outA)
}

plist <- list(
    dat      = list(dfA),
    Y        = 'logitfluP',
    xlag     = 0:2
  ) %>% cross_df()
  
felm_results_all <- plist %>% pmap_df(felm_function_all) %>% 
  as.tibble() %>% 
  mutate(betalow=beta-3.29*se, 
         betahigh=beta+3.29*se) %>% # 99.9% CI
  mutate(state=as.factor(state),
         lag=factor(lag, levels = c(0,1,2),
                    labels = c("Lag 0","Lag 1","Lag 2")),
         sig=ifelse(p<0.001, 'sig', 'non_sig'),
         sig=factor(sig, levels=c("non_sig","sig"))) %>%
  select(state, plt, lag, beta, betalow, betahigh, p, sig)

abc_levels=str_sort(toupper(unique(felm_results_state$state)),decreasing=F)

felm_results <- rbind(felm_results_all, felm_results_state) %>% 
  mutate(state=factor(toupper(state),
                      levels=c(abc_levels,"ALL")))


# Calculate overall and state-specific SD of each environmental predictor
df_SD_all=dfA %>%
  dplyr::summarize(o3=sd(o3), ah=sd(ah), temp=sd(temp)) %>%
  gather(plt, SD) %>% 
  mutate(state='ALL')

df_SD_ST <- dfA %>% group_by(state) %>%
  dplyr::summarize(o3=sd(o3), ah=sd(ah), temp=sd(temp)) %>%
  gather(plt, SD, -state) %>%
  mutate(state=toupper(state))

df_SD <- rbind(df_SD_ST, df_SD_all%>%select(state,plt,SD))

# Transform beta to effect size of each SD change
df_felm_final=felm_results %>%
  left_join(df_SD,by=c("state",'plt')) %>%
  mutate(Size=beta*SD,SizeL=betalow*SD,SizeH=betahigh*SD) %>%
  select(state, plt, lag, Size, SizeL, SizeH, gam.p=p, sig) %>%
  mutate(plt=factor(plt, levels=c( "o3","ah", "temp"),
                    labels=c("O[3]", "AH", "T"))) %>% 
  print(n=Inf)


########################################## Plot ##########################################
mytheme <- theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = NA, colour ="grey90" ),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text.x = element_text(size = 14, color = "gray10",
                                    face='bold',family='serif'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid",
                                 colour = "black"))+
  theme(legend.position = c(.85, 0.025), legend.box = "horizontal",
        legend.title=element_text(size=18,family='serif'),
        legend.key.width= unit(1.1, 'cm'),
        legend.text = element_text(size=14, color = "black",family='serif'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.background = element_rect(color = NA),
        legend.box.margin = margin(0.1,0.1,0.1,0.1,"cm")) +
  theme(axis.title = element_text(size=20,family='serif'),
        axis.text= element_text(color="black", size=18,family='serif'),
        plot.title = element_text(size = 20, hjust=0.5, family='serif'))

# Plotting tate-specific felm results
dat_felm_ST_lag1 <- df_felm_final %>% 
  filter(state!='ALL' & lag=='Lag 1') %>% 
  mutate(state=factor(state, levels= abc_levels)) 

p_felm_ST_lag1= ggplot(dat_felm_ST_lag1) + 
  geom_errorbar(aes(ymin=SizeL, ymax=SizeH, x=fct_rev(plt)), color='gray',
                position = position_dodge(0.8), width=0, size=0.8) +
  geom_point(aes(y=Size, x=fct_rev(plt), shape=sig),
             color='red', size=4, stroke = 0.5,
             position = position_dodge(0.8))+
  facet_wrap(~ state, ncol = 6)+
  scale_shape_manual(name="Statistical significance test:", values = c(1,16),
                     labels=c("Non-significant","Significant")) +
  labs(x='',y=expression(paste(beta, ' estimates in FELM'))) +
  scale_x_discrete(labels = c(expression("T"), expression("AH"), expression(O[3]))) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") +
  coord_flip() + 
  mytheme

# Plotting nationwide FELM results
dat_felm_ALL_lag1 <- df_felm_final %>% 
  filter(state=='ALL' & lag=='Lag 1')

p_felm_vs_lag1 <- ggplot(dat_felm_ALL_lag1) +
  geom_errorbar(aes(ymin=SizeL, ymax=SizeH, x=fct_rev(plt)), color='gray',
                position = position_dodge(0.8), width=0,linewidth=1) +
  geom_point(aes(y=Size, x=fct_rev(plt), shape=sig),
             color='red', size=4, stroke = 0.5,
             position = position_dodge(0.8))+
  scale_shape_manual(guide="none",values = c(1,16)) +
  labs(x='',y=expression(paste(beta, ' estimates in FELM'))) + 
  scale_x_discrete(labels = c(expression("T"), expression("AH"), expression(O[3]))) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") +
  coord_flip() +
  mytheme +
  theme(panel.background = element_blank())

p_felm_vs_lag1_arrange = ggarrange(NULL,p_felm_vs_lag1,NULL,
                               nrow = 1, ncol = 3,
                               widths = c(0.5, 1, 0.5),
                               common.legend = FALSE, 
                               align = "hv")

p_felm_3_all <- ggarrange(p_felm_ST_lag1,
                         p_felm_vs_lag1_arrange,
                         nrow = 2, ncol = 1,
                         heights = c(11, 2.5),
                         common.legend = FALSE,
                         align = "hv")

print(p_felm_3_all)
