################################## Load packages ###################################
options(warn = -1)
options(width=80) 
options(stringsAsFactors=FALSE)
options(scipen = 6) # bias against scientific notation
options(digits = 3) # show fewer decimal places

library(tidyverse); library(lubridate); library(dlnm); library(splines);
library(tsModel); library(gnm);library(ggpubr);library(metafor);library(gam)

################################## Load data ###################################
load('usa_Flu_J_proxy_Data.rda')
dfA=usa_Flu_J_proxy_Data
dfA <- dfA %>% mutate(ym=paste(year(date),month(date),sep = ""))
dfA <- dfA %>% mutate(state=case_when(state == "al" ~ "Alabama", state == "ar" ~ "Arkansas",
                                      state == "az" ~ "Arizona", state == "ca" ~ "California",
                                      state == "co" ~ "Colorado", state == "ct" ~ "Connecticut",
                                      state == "de" ~ "Delaware", state == "ga" ~ "Georgia",
                                      state == "il" ~ "Illinois", state == "ky" ~ "Kentucky",
                                      state == "la" ~ "Louisiana", state == "md" ~ "Maryland",
                                      state == "mn" ~ "Minnesota", state == "mo" ~ "Missouri",
                                      state == "mt" ~ "Montana", state == "ny" ~ "New York",
                                      state == "oh" ~ "Ohio", state == "ok" ~ "Oklahoma",
                                      state == "or" ~ "Oregon", state == "pa" ~ "Pennsylvania",
                                      state == "sc" ~ "South Carolina", state == "sd" ~ "South Dakota",
                                      state == "tx" ~ "Texas", state == "ut" ~ "Utah", 
                                      state == "va" ~ "Virginia", state == "wa" ~ "Washington",
                                      state == "wi" ~ "Wisconsin", state == "wv" ~ "West Virginia",
                                      state == "hi" ~ "Hawaii", state == "ak" ~ "Alaska"))

################################## State-specific GAM analysis ###################################
outdlnm <- function(cityname,lag){
  city <- dfA %>% filter(state == cityname)
  city <- city %>% mutate(lag1 = Lag(fluP, 1)) %>% na.omit()
  # crossbasis of O3
  cb.o3 <<- crossbasis(city$o3, lag = lag, argvar = list(fun = "lin"),
                       arglag = list(fun = 'integer'))
  # model
  model <- gam(fluP ~ cb.o3 + ns(temp,6) + ns(lag(temp,1),6) +
                 ns(ah,6) + ns(lag(ah,1),6) + lag1 +  as.factor(ym) ,
               family = quasibinomial, data =city)
  iqr <- IQR(city$temp, na.rm=T)
  pred <- crosspred(cb.o3, model, at = iqr, cumul = T)
  beta <- rownames_to_column(data.frame(pred$allfit), var = "conc")
  se <- rownames_to_column(data.frame(pred$allse), var = "conc")
  out <- left_join(beta, se, by = "conc") %>% mutate(conc = as.numeric(conc))
  out$city <- cityname
  out$lag <- lag
  out
}

plist=list(cityname=unique(dfA$state),lag=0:2) %>% cross_df()
out <- plist %>% pmap_df(outdlnm)
names(out) <- c("conc", "allfit", "allse",  "state", "lag")
output <- out %>% select(state,allfit,allse,lag) %>% mutate(allsel=allfit-1.96*allse,allseh=allfit+1.6*allse)
output$p <-  pnorm(-abs(out$allfit) / out$allse)  * 2
GAMres <- output %>% select(state,lag,allfit,p) %>% mutate(plt="o3")
names(GAMres) <- c("state","lag","effect","P_value","plt")

################################## Meta analysis of state-specific results ###################################
out1 <- out %>% filter(lag==0)
meta.re <- rma(yi=allfit, sei=allse, slab=state, method="REML", data=out1)
result1 <- with(meta.re, c(b, se, pval))

out1 <- out %>% filter(lag==1)
meta.re <- rma(yi=allfit, sei=allse, slab=state, method="REML", data=out1)
result2 <- with(meta.re, c(b,se, pval))

out1 <- out %>% filter(lag==2)
meta.re <- rma(yi=allfit, sei=allse, slab=state, method="REML", data=out1)
result3 <- with(meta.re, c(b, se, pval))

result <- rbind(result1,result2,result3) %>% as.data.frame() %>% mutate(state="OVERALL",lag=c(0:2),p=V3)
result <- result %>% select(state,lag,V1,p) %>% mutate(plt="o3")
names(result) <- c("state","lag","effect","P_value","plt")

################################## All the results: state-specific + overall meta ###################################
GAMres <- rbind(GAMres,result) 

################################## plot ###################################
output <- out %>% select(state,allfit,allse,lag) %>% mutate(allsel=allfit-allse*1.96,allseh=allfit+allse*1.96) 
output <- output[order(output$allfit),]
output$lag <- factor(output$lag, levels = 0:2, labels = c("Lag 0","Lag 1","Lag 2"))
output$state <- as.factor(output$state)
output <- output %>% select("state","lag","allfit","allsel","allseh")

result <- as.data.frame(matrix(rep(NA,15),nrow=3))
names(result) <- names(output)

out1 <- out %>% filter(lag==0)
meta.re <- rma(yi=allfit, sei=allse, slab=state, method="REML", data=out1)
result1 <- with(meta.re, c(b, ci.lb, ci.ub))
a <- result1[1]
b <- result1[2]
c <- result1[3]
result[1,3:5] <- c(a,b,c)

out1 <- out %>% filter(lag==1)
meta.re <- rma(yi=allfit, sei=allse, slab=state, method="REML", data=out1)
result2 <- with(meta.re, c(b, ci.lb, ci.ub))
a <- result2[1]
b <- result2[2]
c <- result2[3]
result[2,3:5] <- c(a,b,c)

out1 <- out %>% filter(lag==2)
meta.re <- rma(yi=allfit, sei=allse, slab=state, method="REML", data=out1)
result3 <- with(meta.re, c(b, ci.lb, ci.ub))
a <- result3[1]
b <- result3[2]
c <- result3[3]
result[3,3:5] <- c(a,b,c)

result$state <- "Overall"
result$lag <- 0:2
result$lag <- factor(result$lag, levels = 0:2, labels = c("Lag 0","Lag 1","Lag 2"))
output <- rbind(output,result) %>% mutate(type=ifelse(state=="Overall",2,1))
output$type <- factor(output$type)

mytheme3 <-  theme_bw() +
  theme(axis.title = element_text(size=18),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(vjust= 2.5),
        axis.text.x=element_text(size = 16, colour = "black"),
        axis.text.y=element_text(size = 16, colour = "black"),
        axis.ticks.length=unit(0.1,'cm'),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position="none",
        strip.text = element_text(size = rel(1.2), face = "bold"),
        plot.title = element_text(size=18, hjust=0.5),
        plot.margin = unit(c(0,0.1,0.5,0.5), "lines"))

ggplot(output) + 
  geom_point(aes(y=allfit, x=state,col=type))+
  geom_linerange(aes(ymin=allsel, ymax=allseh, x=state,col=type)) +
  scale_color_manual(values = c("black","red")) +
  geom_hline(yintercept = 1, col="red",linetype = 2) +
  coord_flip() +
  xlab("") +
  ylab("Causal effect estimates by GAM") +
  facet_grid(~lag,scales='free')+
  mytheme3
