################################## Load packages & data ###################################
library(tidyverse); library(lubridate); library(dlnm); library(splines);
library(tsModel); library(gnm);library(ggpubr);library(metafor);
library(mgcv); library(imputeTS)

load("Data_fluseason.rda")
df <- usa_Flu_P_proxy_Data_B

dfA <- df %>% 
  mutate(time = as.numeric(date), 
         year = year(date),
         month = month(date))

################################## State-level GLM ###################################
gam_function <- function(dat, Y='fluP', x1='value', xlag=0, cityname, variable) {
  
  df <- dat %>%
    filter(state == cityname) %>%
    mutate(loglag1 = log(Lag(fluP, 1)),
           loglag2 = Lag(loglag1, 1))
  
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

########################### Meta analysis of GLM beta results ############################

fmeta=function(dat){
  out=dat
  
  GAMresCI  =  out %>%
    mutate(betalow=beta-1.96*se, 
           betahigh=beta+1.96*se) %>% # 95% CI
    mutate(state=as.factor(state),
           lag=factor(lag, levels = c(0,1,2),
                      labels = c("Lag 0","Lag 1","Lag 2"))) %>%
    select(state,lag,beta,betalow,betahigh,p)
  
  meta.ci.result  =  as.data.frame(matrix(rep(NA,18),nrow=3))
  names(meta.ci.result)  =  names(GAMresCI)
  
  for (i in 1:3) {
    out.lag  =  out %>%
      filter(lag==i-1)
    meta  =  rma(yi=beta, sei=se, slab=state, method="REML", data=out.lag, level=99.9)
    meta.re  =  with(meta, c(b, ci.lb, ci.ub,pval))
    meta.ci.result[i,1] = "All"
    meta.ci.result[i,2] = paste('Lag',i-1)
    meta.ci.result[i,3] = meta.re[1]
    meta.ci.result[i,4] = meta.re[2]
    meta.ci.result[i,5] = meta.re[3]
    meta.ci.result[i,6] = meta.re[4]
  }
  
  meta.ci.result  =  meta.ci.result %>% as.data.frame() 
  
  abc_levels=str_sort(toupper(unique(GAMresCI$state)),decreasing=T)
  
  GAMresCI_all  =  rbind(GAMresCI, meta.ci.result) %>%
    mutate(type=factor(ifelse(state=="All",2,1))) %>%
    mutate(state=factor(toupper(state),
                        levels=c(abc_levels,"ALL"))) %>%
    as.tibble()
}

gam_o3=fmeta(out_o3)
gam_ah=fmeta(out_ah)
gam_temp=fmeta(out_temp)

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
df1=gam_o3  %>% as.data.frame() %>% mutate(plt='o3')
df2=gam_ah  %>% as.data.frame() %>% mutate(plt='ah')
df3=gam_temp %>% as.data.frame() %>% mutate(plt='temp')

df_gam_final=rbind(df1,df2,df3) %>%
  left_join(df_SD,by=c("state",'plt')) %>%
  mutate(Size=beta*SD,SizeL=betalow*SD,SizeH=betahigh*SD) %>%
  select(state, plt, lag, Size, SizeL, SizeH, gam.p=p) %>%
  mutate(sig=ifelse(gam.p<0.05, 'sig', 'non_sig'),
         sig=factor(sig, levels=c("non_sig","sig")),
         plt=factor(plt, levels=c('temp', 'ah', 'o3')))

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

# Plotting tate-specific GLM results
abc_levels=str_sort(unique(df_gam_final$state)[1:46],decreasing=F)

dat_glm_ST <- df_gam_final %>% 
  filter(state!='ALL' & lag=='Lag 1') %>% 
  mutate(state=factor(state, levels= abc_levels)) 

p_glm_ST_lag1= ggplot(dat_glm_ST) + 
  geom_errorbar(aes(ymin=SizeL, ymax=SizeH, x=plt),
                position = position_dodge(0.8), width=0, size=0.8) +
  geom_point(aes(y=Size, x=plt, shape=sig),
             fill='black', size=4, stroke = 0.5,
             position = position_dodge(0.8))+
  facet_wrap(~ state, ncol = 6)+
  scale_shape_manual(name="Statistical significance test:", values = c(1,16),
                     labels=c("Non-significant","Significant")) +
  labs(x='',y=expression(paste(beta, ' estimates in GLM'))) +
  scale_x_discrete(labels = c(expression("T"), expression("AH"), expression(O[3]))) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") +
  coord_flip() +
  mytheme

# Plotting meta-analyzed GLM results
dat_glm_ALL <- df_gam_final %>% 
  filter(state=='ALL' & lag=='Lag 1')

p_glm_vs_1 <- ggplot(dat_glm_ALL) +
  geom_errorbar(aes(ymin=SizeL, ymax=SizeH, x=plt),
                position = position_dodge(0.8), width=0,linewidth=1) +
  geom_point(aes(y=Size, x=plt, shape=sig),
             fill='black', size=4, stroke = 0.5,
             position = position_dodge(0.8))+
  scale_shape_manual(guide="none",values = c(1,16)) +
  labs(x='',y=expression(paste(beta, ' estimates in GLM'))) + 
  scale_x_discrete(labels = c(expression("T"), expression("AH"), expression(O[3]))) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") +
  coord_flip() +
  mytheme

p_glm_vs_1_arrange = ggarrange(NULL,p_glm_vs_1,NULL,
                               nrow = 1, ncol = 3,
                               widths = c(0.5, 1, 0.5),
                               common.legend = FALSE, 
                               align = "hv")

p_glm_3_all <- ggarrange(p_glm_ST_lag1,
                         p_glm_vs_1_arrange,
                         nrow = 2, ncol = 1,
                         heights = c(11, 2.5),
                         common.legend = FALSE,
                         align = "hv")

#+ fig.height=12,fig.width=10
p_glm_3_all 
#+
