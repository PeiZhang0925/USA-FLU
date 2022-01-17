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
load("usa_Flu_M_proxy_Data.rda")

dfA=usa_Flu_M_proxy_Data
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
  rename(ST=state)

dfA <- dfA %>% mutate(state=case_when(ST == "al" ~ "Alabama", ST == "ar" ~ "Arkansas",
                                   ST == "az" ~ "Arizona", ST == "ca" ~ "California",
                                   ST == "co" ~ "Colorado", ST == "ct" ~ "Connecticut",
                                   ST == "de" ~ "Delaware", ST == "ga" ~ "Georgia",
                                   ST == "il" ~ "Illinois", ST == "ky" ~ "Kentucky",
                                   ST == "la" ~ "Louisiana", ST == "md" ~ "Maryland",
                                   ST == "mn" ~ "Minnesota", ST == "mo" ~ "Missouri",
                                   ST == "mt" ~ "Montana", ST == "ny" ~ "New York",
                                   ST == "oh" ~ "Ohio", ST == "ok" ~ "Oklahoma",
                                   ST == "or" ~ "Oregon", ST == "pa" ~ "Pennsylvania",
                                   ST == "sc" ~ "South Carolina", ST == "sd" ~ "South Dakota",
                                   ST == "tx" ~ "Texas", ST == "ut" ~ "Utah", 
                                   ST == "va" ~ "Virginia", ST == "wa" ~ "Washington",
                                   ST == "wi" ~ "Wisconsin", ST == "wv" ~ "West Virginia",
                                   ST == "hi" ~ "Hawaii", ST == "ak" ~ "Alaska",
                                   ST == "tn" ~'Tennessee'))

################################## State-specific GAM analysis ###################################

fgam2=function(Y='cvsd',x1='value',x2,x3,xlag=0,cityname){
  
  df <- dfA %>% filter(state == cityname) %>% 
    mutate(loglag1 = log(Lag(fluP, 1)+0.0000001)) 
  
  fit=gam(as.formula(paste(Y,
                           '~Lag(',x1,',',xlag,')+Lag(',x2,',',xlag,')+
                           Lag(',x3,',',xlag,')+
                           as.factor(year) + as.factor(month)+loglag1',sep='')),
          data=df,family=quasibinomial,na.action=na.omit, 
          #weights=w,
          control = glm.control(epsilon = 1e-10, maxit = 1000))
  
  summ=summary(fit)
  
  outA=data.frame(summ$p.coeff[2],summ$se[2],summ$p.t[2],summ$p.pv[2])
  names(outA)=c('beta','se','t','p')
  outA$state=cityname
  outA$lag=xlag
  outA
}

plist=list(Y='fluP',x1='o3',x2='temp',x3='ah',
           xlag=0:2,
           cityname=unique(dfA$state)) %>% cross_df()
out <- plist %>% pmap_df(fgam2)

out

GAMres <- out %>% select(state,lag, effect=beta, p) %>% as.tibble()


################################## Meta analysis of state-specific results ###################################

result <- matrix(nrow=3, ncol=6, dimnames=list(1:3, c("state", "lag","effect",'p','lb','ub')))

for (i in 1:3) {
  out.lag <- out %>% filter(!state %in% c("Hawaii","Alaska"))%>% filter(lag==i-1)
  meta <- rma(yi=beta, sei=se, slab=state, method="REML", data=out.lag)
  meta.re <- with(meta, c(b, pval, ci.lb,ci.ub))
  result[i,1] = "OVERALL"
  result[i,2] = i-1
  result[i,3] = meta.re[1]
  result[i,4] = meta.re[2]
  result[i,5] = meta.re[3]
  result[i,6] = meta.re[4]
}

meta.result <- result %>%
  as.data.frame() %>%
  mutate_at(3:6, as.numeric)

meta.result

################################## All the results: state-specific + overall meta ###################################

roundFunc <- function(x){round(x, 3)}
GAMres <-GAMres %>% filter(!state %in% c("Hawaii","Alaska"))
order0 <- GAMres[GAMres$lag==0,]
order0 <- order0[order(order0$effect),]
order1 <- GAMres[GAMres$lag==1,]
order1 <- order1[order(order1$effect),]
order2 <- GAMres[GAMres$lag==2,]
order2 <- order2[order(order2$effect),]


################################## plot ###################################

head(out) 

GAMresCI <- out %>% 
  mutate(betalow=beta-1.96*se,
         betahigh=beta+1.96*se) %>%
  mutate(state=as.factor(state),
         lag=factor(lag, levels = c(0,1,2),
                    labels = c("Lag 0","Lag 1","Lag 2"))) %>%
  select(state,lag,beta,betalow,betahigh)

ci.result <- as.data.frame(matrix(rep(NA,15),nrow=3))
names(ci.result) <- names(GAMresCI)
# state lag beta betalow betahigh

for (i in 1:3) {
  out.lag <- out %>% filter(!state %in% c("Hawaii","Alaska"))%>% filter(lag==i-1)
  meta <- rma(yi=beta, sei=se, slab=state, method="REML", data=out.lag)
  meta.re <- with(meta, c(b, ci.lb, ci.ub))
  ci.result[i,1] = "Overall"
  ci.result[i,2] = paste('Lag',i-1)
  ci.result[i,3] = meta.re[1]
  ci.result[i,4] = meta.re[2]
  ci.result[i,5] = meta.re[3]
}

GAMresCI_all <- rbind(GAMresCI,ci.result) %>% 
  mutate(type=factor(ifelse(state=="Overall",2,1))) %>%
  as.tibble()

mytheme3 <-  theme_bw() +
  theme(axis.title = element_text(size=14),
        axis.title.x = element_text(margin = margin(t = 18)),
        axis.title.y = element_text(vjust= 2.5),
        axis.text.x=element_text(size = 12, colour = "black"),
        axis.text.y=element_text(size = 12, colour = "black"),
        axis.ticks.length=unit(0.1,'cm'),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position="none",
        strip.text = element_text(size = rel(1.2), face = "bold"),
        plot.title = element_text(size=18, hjust=0.5),
        plot.margin = unit(c(0.3,0.3,0,0), "lines"))

mytheme2 <- theme_bw() + 
  theme(axis.text.x=element_text(size = 12, color= "black"),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.ticks.length=unit(0.1,'cm'),
        panel.grid = element_blank(),
        legend.position="none",
        plot.title = element_text(size = 18, hjust=0.5))
  
GAMresCI_all <- GAMresCI_all %>%  filter(!state %in% c("Hawaii","Alaska"))

dat0 <- GAMresCI_all[GAMresCI_all$lag=="Lag 0",]

dat0$state <- factor(dat0$state,levels = c(order1$state,"Overall"),
                      labels = c(order1$state,"Overall"))

p1 <- ggplot(dat0) + 
  geom_point(aes(y=beta, x=state,col=type))+
  # geom_linerange(aes(ymin=betalow, ymax=betahigh, x=state,col=type)) +
  geom_errorbar(aes(ymin=betalow, ymax=betahigh, x=state,col=type, width=0.3)) +
  scale_color_manual(values = c("black","blue")) +
  scale_x_discrete(labels=c(c(order1$state,expression(bold("Overall"))))) +
  geom_hline(yintercept = 0, col="red",linetype = 2) +
  coord_flip() +
  ggtitle("Lag 0") +
  labs(x='',y='') +
  mytheme3


dat1 <- GAMresCI_all[GAMresCI_all$lag=="Lag 1",]

dat1$state <- factor(dat1$state,levels = c(order1$state,"Overall"),
                     labels =  c(order1$state,"Overall"))

p2 <- ggplot(dat1) + 
  geom_point(aes(y=beta, x=state,col=type))+
  geom_errorbar(aes(ymin=betalow, ymax=betahigh, x=state,col=type, width=0.3)) +
  scale_color_manual(values = c("black","blue")) +
  scale_x_discrete(labels=c(c(order1$state,expression(bold("Overall"))))) +
  geom_hline(yintercept = 0, col="red",linetype = 2) +
  coord_flip() +
  ggtitle("Lag 1") +
  labs(x='',y='') +
  mytheme2

dat2 <- GAMresCI_all[GAMresCI_all$lag=="Lag 2",]

dat2$state <- factor(dat2$state,levels = c(order1$state,"Overall"),
                     labels =  c(order1$state,"Overall"))

p3 <- ggplot(dat2) + 
  geom_point(aes(y=beta, x=state,col=type))+
  geom_errorbar(aes(ymin=betalow, ymax=betahigh, x=state,col=type, width=0.3)) +
  scale_color_manual(values = c("black","blue")) +
  scale_x_discrete(labels=c(c(order1$state,expression(bold("Overall"))))) +
  geom_hline(yintercept = 0, col="red",linetype = 2) +
  coord_flip() +
  ggtitle("Lag 2") +
  labs(x='',y='') +
  mytheme2

fig2 <- ggarrange(p1, p2, p3,
                  nrow = 1, ncol = 3, widths = c(1.47,1,1),align = "h") %>%
  annotate_figure(bottom  = text_grob(expression(paste("Change in Log Odds (Flu) per ppb ", O[3], " increse")), 
                                      size = 14, face="bold", vjust =-1.8, hjust=0.4))

#+ fig.height=10,fig.width=9, out.width = "120%"
fig2 
#+
ggsave(fig2,file="GAM_final_newFunc_yrFct.jpg", width = 8, height = 8, dpi = 500)


########################sensitive analysis######################
fgam2=function(Y='cvsd',x1='value',x2,x3,xlag=0,cityname){
  
  df <- dfA %>% filter(state == cityname) %>% 
    mutate(loglag1 = log(Lag(fluP, 1)+0.0000001)) 
  
  fit=gam(as.formula(paste(Y,
                           # '~Lag(',x1,',',xlag,')+ns(', x2,',6)+ns(Lag(',x2,',1),6)+ ns(',x3,',6)+ns(Lag(',x3,',1),6)+ as.factor(year) + as.factor(month)+loglag1',
                           '~Lag(',x1,',',xlag,')+ns(Lag(',x2,',', xlag,'),6)+ ns(Lag(',x3,',', xlag,'),6)+ s(time,k=(5*8+1),fx=T)+loglag1',
                           sep='')),
          data=df,family=quasibinomial,na.action=na.omit, 
          #weights=w, control = glm.control(epsilon = 1e-10, maxit = 1000)
  )
  
  summ=summary(fit)
  
  outA=data.frame(summ$p.coeff[2],summ$se[2],summ$p.t[2],summ$p.pv[2])
  names(outA)=c('beta','se','t','p')
  outA$state=cityname
  outA$lag=xlag
  outA
}

plist=list(Y='fluP',x1='o3',x2='temp',x3='ah',
           xlag=0:2,
           cityname=unique(dfA$state)) %>% cross_df()
out <- plist %>% pmap_df(fgam2)

out

GAMres <- out %>% select(state,lag, effect=beta, p) %>% as.tibble()


################################## Meta analysis of state-specific results ###################################

result <- matrix(nrow=3, ncol=6, dimnames=list(1:3, c("state", "lag","effect",'p','lb','ub')))

for (i in 1:3) {
  out.lag <- out %>% filter(!state %in% c("Hawaii","Alaska"))%>% filter(lag==i-1)
  meta <- rma(yi=beta, sei=se, slab=state, method="REML", data=out.lag)
  meta.re <- with(meta, c(b, pval, ci.lb,ci.ub))
  result[i,1] = "OVERALL"
  result[i,2] = i-1
  result[i,3] = meta.re[1]
  result[i,4] = meta.re[2]
  result[i,5] = meta.re[3]
  result[i,6] = meta.re[4]
}

meta.result <- result %>%
  as.data.frame() %>%
  mutate_at(3:6, as.numeric)

meta.result

################################## All the results: state-specific + overall meta ###################################

GAMres <-GAMres %>% filter(!state %in% c("Hawaii","Alaska"))
order0 <- GAMres[GAMres$lag==0,]
order0 <- order0[order(order0$effect),]
order1 <- GAMres[GAMres$lag==1,]
order1 <- order1[order(order1$effect),]
order2 <- GAMres[GAMres$lag==2,]
order2 <- order2[order(order2$effect),]


################################## plot ###################################

head(out) 

GAMresCI <- out %>% 
  mutate(betalow=beta-1.96*se,
         betahigh=beta+1.96*se) %>%
  mutate(state=as.factor(state),
         lag=factor(lag, levels = c(0,1,2),
                    labels = c("Lag 0","Lag 1","Lag 2"))) %>%
  select(state,lag,beta,betalow,betahigh)

ci.result <- as.data.frame(matrix(rep(NA,15),nrow=3))
names(ci.result) <- names(GAMresCI)
# state lag beta betalow betahigh

for (i in 1:3) {
  out.lag <- out %>% filter(!state %in% c("Hawaii","Alaska"))%>% filter(lag==i-1)
  meta <- rma(yi=beta, sei=se, slab=state, method="REML", data=out.lag)
  meta.re <- with(meta, c(b, ci.lb, ci.ub))
  ci.result[i,1] = "Overall"
  ci.result[i,2] = paste('Lag',i-1)
  ci.result[i,3] = meta.re[1]
  ci.result[i,4] = meta.re[2]
  ci.result[i,5] = meta.re[3]
}

GAMresCI_all <- rbind(GAMresCI,ci.result) %>% 
  mutate(type=factor(ifelse(state=="Overall",2,1))) %>%
  as.tibble()

GAMresCI_all <- GAMresCI_all %>%  filter(!state %in% c("Hawaii","Alaska"))

dat0 <- GAMresCI_all[GAMresCI_all$lag=="Lag 0",]

dat0$state <- factor(dat0$state,levels = c(order1$state,"Overall"),
                     labels = c(order1$state,"Overall"))

p1 <- ggplot(dat0) + 
  geom_point(aes(y=beta, x=state,col=type))+
  # geom_linerange(aes(ymin=betalow, ymax=betahigh, x=state,col=type)) +
  geom_errorbar(aes(ymin=betalow, ymax=betahigh, x=state,col=type, width=0.3)) +
  scale_color_manual(values = c("black","blue")) +
  scale_x_discrete(labels=c(c(order1$state,expression(bold("Overall"))))) +
  geom_hline(yintercept = 0, col="red",linetype = 2) +
  coord_flip() +
  ggtitle("Lag 0") +
  labs(x='',y='') +
  mytheme3


dat1 <- GAMresCI_all[GAMresCI_all$lag=="Lag 1",]

dat1$state <- factor(dat1$state,levels = c(order1$state,"Overall"),
                     labels =  c(order1$state,"Overall"))

p2 <- ggplot(dat1) + 
  geom_point(aes(y=beta, x=state,col=type))+
  geom_errorbar(aes(ymin=betalow, ymax=betahigh, x=state,col=type, width=0.3)) +
  scale_color_manual(values = c("black","blue")) +
  scale_x_discrete(labels=c(c(order1$state,expression(bold("Overall"))))) +
  geom_hline(yintercept = 0, col="red",linetype = 2) +
  coord_flip() +
  ggtitle("Lag 1") +
  labs(x='',y='') +
  mytheme2

dat2 <- GAMresCI_all[GAMresCI_all$lag=="Lag 2",]

dat2$state <- factor(dat2$state,levels = c(order1$state,"Overall"),
                     labels =  c(order1$state,"Overall"))

p3 <- ggplot(dat2) + 
  geom_point(aes(y=beta, x=state,col=type))+
  geom_errorbar(aes(ymin=betalow, ymax=betahigh, x=state,col=type, width=0.3)) +
  scale_color_manual(values = c("black","blue")) +
  scale_x_discrete(labels=c(c(order1$state,expression(bold("Overall"))))) +
  geom_hline(yintercept = 0, col="red",linetype = 2) +
  coord_flip() +
  ggtitle("Lag 2") +
  labs(x='',y='') +
  mytheme2

fig2 <- ggarrange(p1, p2, p3,
                  nrow = 1, ncol = 3, widths = c(1.47,1,1),align = "h") %>%
  annotate_figure(bottom  = text_grob(expression(paste("Change in Log Odds (Flu) per ppb ", O[3], " increse")), 
                                      size = 14, face="bold", vjust =-1.8, hjust=0.4))

#+ fig.height=10,fig.width=9, out.width = "120%"
fig2 
#+
ggsave(fig2,file="GAM_final_newFunc_yrFct_sensitive.jpg", width = 8, height = 8, dpi = 500)
