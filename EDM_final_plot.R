
################################## Load packages ###################################

# opts_chunk$set(warning=FALSE,message=FALSE)
options(warn = -1)
options(width=80) 
options(stringsAsFactors=FALSE)
options(scipen = 6) # bias against scientific notation
options(digits = 3) # show fewer decimal places

# packages
require(lubridate); require(tidyverse)
require(ggplot2);require(ggpubr)
require(ggthemes);require(cowplot)
require(customLayout); require(patchwork)
require(grid); require(gridExtra)
require(usmap);require(maps)
require(quantreg)
require(metap)

OutputPath <- paste0("E:/Dropbox (linweit)/ehrg/users/Fang/US_FLU/plot_NCappeal_reproduce_H/")

################################## Load data ###################################

load('usa_Flu_J_proxy_Data.rda')

dfA=usa_Flu_J_proxy_Data

normFunc=function(x){(x-mean(x, na.rm = T))/sd(x, na.rm = T)}

df_smapc=dfA %>% group_by(state) %>% 
  mutate_at(3:ncol(.),normFunc) %>% ungroup()

df_flu=dfA %>% select(date,state,ah,o3,temp,fluR,fluP) %>% group_by(state) %>% 
  gather(key,value,-c(date,state))

df_T <- dfA %>%
  select(date, state, temp, ah, o3, fluR, fluP)

ts <- df_T %>%
  select(-fluR) %>%
  gather(key,value,-date,-state) %>%
  mutate(variable= factor(key,levels = c("ah", "temp", "o3", "fluP"),
                          labels = c("AH", "T", "O[3]", "Influenza"))) %>%
  ggplot(aes(date, value, color = variable)) +
  geom_point(size=1) +
  facet_grid(variable ~ ., scales = 'free', labeller = label_parsed) +
  scale_x_date(breaks = seq(as.Date("2011-01-01"),
                            as.Date("2015-01-01"), by="1 year"),
               date_labels="%b %Y") +
  labs(x="Week", y="Values") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 12),
    axis.title.y = element_text(vjust= 2.5),
    axis.title.x = element_text(margin = margin(t = 12)),
    axis.text.x=element_text(size = 10, color ="black"),
    axis.text.y=element_text(size = 10, color ="black"),
    axis.ticks.length=unit(0.1,'cm'),
    legend.position="none",
    panel.grid = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(size = rel(1.3), face = "bold")
  )

ts

#################################### 1. CCM_OUT: perform CCM for real and surrogate data ###################################

##### 1) Calculate p-value for CCM test #####

load('ccm_out_2021flu_o3_H_alpha.rda')

unique(ccm_out$ST)

ccm_p=ccm_out %>% 
  group_by(dis,plt,E,ST,tp_value) %>%
  summarise(p=1-ecdf(rho[i != 1])(rho[i == 1]))

print(ccm_p, n=500)

fn_metap_all=function(tp_values,plts,diss){
  df=filter(ccm_p,tp_value==tp_values & plt==plts & dis==diss)
  out=allmetap(df$p, method = "all") %>% as.data.frame()
  mutate(out,tp_value=tp_values,plt=plts,dis=diss)
}

plist=list(tp_values=-2:0, 
           plts=c('o3',"temp","ah"),
           diss='fluP') %>% cross_df()

ccm_metap_all=plist %>% pmap_df(fn_metap_all)

ccm_metap_all

fn_metap=function(tp_values,plts,diss){
  df=filter(ccm_p,tp_value==tp_values & plt==plts & dis==diss)
  out=allmetap(df$p, method = "sumlog") %>% as.data.frame()
  mutate(out,tp_value=tp_values,plt=plts,dis=diss)
}

metap_out=plist %>% pmap_df(fn_metap)

metap_out

##### 2) Plotting CCM test #####

names(ccm_out); names(ccm_p)
sort(unique(ccm_out$ST)); sort(unique(ccm_p$ST))

ccm_causal <- ccm_out %>%
  full_join(ccm_p, by = c("ST","plt","dis",'E','tp_value')) %>%
  mutate(ST=case_when(ST == "al" ~ "Alabama", ST == "ar" ~ "Arkansas",
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
                      ST == "hi" ~ "Hawaii", ST == "ak" ~ "Alaska")) %>%
  drop_na()
names(ccm_causal); unique(ccm_causal$ST); unique(ccm_causal$plt)


##############

mytheme1 <- theme_bw() +
  theme(axis.text.x=element_text(size = 16, color= "black"),
        axis.text.y=element_text(size = 16, color= "black"),
        axis.ticks.length=unit(0.1,'cm'),
        legend.position="none",
        plot.title = element_text(size = 18, hjust=0.5))

mytheme2 <- theme_bw() + 
  theme(axis.text.x=element_text(size = 16, color= "black"),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.ticks.length=unit(0.1,'cm'),
        legend.position="none",
        plot.title = element_text(size = 18, hjust=0.5))

fn_plt_ccm <- function(lag){
  
  p1_fluP <- ccm_causal %>% filter(tp_value==lag) %>%
    filter(dis=="fluP", plt=="ah") %>% 
    mutate(grp=ifelse(i==1,'raw','surr'), sig=as.factor(ifelse(p<=0.05, 'sig', 'non_sig'))) %>%
    mutate(ST=factor(ST, levels = c(rev(unique(ST))))) %>%
    spread(grp,rho)%>% 
    ggplot() + 
    geom_boxplot(aes(ST,surr)) +
    geom_point(aes(x=ST,y=raw, shape=sig),size=2, color='red') +
    scale_shape_manual(values = c(1, 16)) +
    labs(x='',y='') +
    ggtitle("AH") + 
    scale_y_continuous(limits = c(-0.325, 0.905), breaks = seq(-0.3, 0.9, by= 0.3)) +
    coord_flip() + mytheme1
  
  p2_fluP <- ccm_causal %>% filter(tp_value==lag) %>%
    filter(dis=="fluP", plt=="temp") %>% 
    mutate(grp=ifelse(i==1,'raw','surr'), sig=as.factor(ifelse(p<=0.05, 'sig', 'non_sig'))) %>%
    mutate(ST=factor(ST, levels = c(rev(unique(ST))))) %>%
    spread(grp,rho)%>% 
    ggplot() + 
    geom_boxplot(aes(ST,surr)) +
    geom_point(aes(x=ST,y=raw, shape=sig),size=2, color='red') +
    scale_shape_manual(values = c(1, 16)) +
    labs(x='',y='') +
    ggtitle("T") + 
    scale_y_continuous(limits = c(-0.325, 0.905), breaks = seq(-0.3, 0.9, by= 0.3)) +
    coord_flip() + mytheme2 
  
  p3_fluP <- ccm_causal %>% filter(tp_value==lag) %>%
    filter(dis=="fluP", plt=="o3")%>% 
    mutate(grp=ifelse(i==1,'raw','surr'), sig=as.factor(ifelse(p<=0.05, 'sig', 'non_sig'))) %>%
    mutate(ST=factor(ST, levels = c(rev(unique(ST))))) %>%
    spread(grp,rho)%>% 
    ggplot() + 
    geom_boxplot(aes(ST,surr)) +
    geom_point(aes(x=ST,y=raw, shape=sig),size=2, color='red') +
    scale_shape_manual(values = c(1, 16)) +
    labs(x='',y='') +
    ggtitle(expression(O[3]))  + 
    scale_y_continuous(limits = c(-0.325, 0.905), breaks = seq(-0.3, 0.9, by= 0.3)) +
    coord_flip() + mytheme2 
  
  p_causal_fluP <- ggarrange(p1_fluP, p2_fluP, p3_fluP, nrow = 1, ncol = 3, 
                             widths = c(1.5, 1, 1), align = "h") %>%
    annotate_figure(bottom  = text_grob(expression(rho["CCM"]), 
                                        size = 18, face="bold"))
  
  ggsave(p_causal_fluP, file=paste0(OutputPath, 
                                    "fluP_causal_", 
                                    as.character(abs(lag)),
                                    ".tiff"), 
         width=14, height=14, units = "in", scale = 0.8)
  
  print(p_causal_fluP) 
}

fn_plt_ccm(0) 

fn_plt_ccm(-1) 

fn_plt_ccm(-2) 

ccm_causal %>% filter(dis=="fluP", plt=="o3") %>% summarise(rho=quantile(rho, c(0,1)), prob =c(0,1))


################################# 2. MFI: Multivariate Forecast improvement ##########################################

##### 1) Calculate p-value for MFI test #####

###### to check whether delta_rho>0 for each state ###### 

load('mfi_surr_out_2021flu_o3_H_alpha.rda')

load('theta_out_2021flu_o3_G.rda')

find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pk <- unlist(pks)
  pks <- ifelse(is.null(pk)==T, which.max(x), pk)
  pks
}

best_theta=theta_out %>%
  group_by(state) %>%
  mutate(row_number=1:n()) %>%
  filter(row_number==find_peaks(rho,m=0)[1]) %>%
  as.data.frame()

dfg2 = best_theta %>%
  select(rho_uni=rho, ST=state) %>%
  right_join(mfi_surr_out, by= 'ST') %>%
  group_by(ST, plt_1, tp_value)%>%
  mutate(E= 10) %>%
  mutate(i=i-1,grp=ifelse(i==0,'raw','surr'),
         rho=rho-rho_uni,
         plt=plt_1) 

dfg2.p=dfg2 %>%
  group_by(dis,plt,E,ST,tp_value) %>%
  summarise(p=1-ecdf(rho[i != 0])(rho[i == 0])) %>%
  arrange(dis,plt)

print(dfg2.p, n=Inf)  # comparing raw and surrogate time series

fn_meta_mfi_p=function(tp_values,plts,diss){
  df=filter(dfg2.p,tp_value==tp_values & plt==plts & dis==diss)
  out=allmetap(df$p, method = "sumlog") %>% as.data.frame
  mutate(out,tp_value=tp_values,plt=plts,dis=diss)
}

plist=list(tp_value=-2:0,
           plt=c("o3","temp","ah"),
           dis=c('fluP')) %>% cross_df()

meta_mfi_p_out=plist %>% pmap_df(fn_metap)

meta_mfi_p_out  # meta-p (surrogate significance test) for each predictor

mfi_causal <- dfg2 %>%
  select(i,rho,ST,dis,plt,E,tp_value,grp) %>%
  full_join(dfg2.p, by = c("ST","plt","dis",'E','tp_value')) %>%
  mutate(ST=case_when(ST == "al" ~ "Alabama", ST == "ar" ~ "Arkansas",
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
                      ST == "hi" ~ "Hawaii", ST == "ak" ~ "Alaska")) %>%
  drop_na()
names(mfi_causal); unique(mfi_causal$ST); unique(mfi_causal$plt)


##############

fn_plt_mfi <- function(lag){
  
  p1_fluP <- mfi_causal %>% filter(tp_value==lag) %>%
    filter(dis=="fluP", plt=="ah") %>% 
    mutate(sig=as.factor(ifelse(p<=0.05, 'sig', 'non_sig'))) %>%
    mutate(ST=factor(ST, levels = str_sort(unique(mfi_causal$ST),decreasing=T)))%>%
    spread(grp,rho)%>% 
    ggplot() + 
    geom_boxplot(aes(ST,surr)) +
    geom_point(aes(x=ST,y=raw, shape=sig),size=2, color='red') +
    scale_shape_manual(values = c(1, 16)) +
    labs(x='',y='') +
    ggtitle("AH") + 
    scale_y_continuous(limits = c(-0.15, 0.1)) +
    coord_flip() + mytheme1
  
  p2_fluP <- mfi_causal %>% filter(tp_value==lag) %>%
    filter(dis=="fluP", plt=="temp") %>% 
    mutate(sig=as.factor(ifelse(p<=0.05, 'sig', 'non_sig'))) %>%
    mutate(ST=factor(ST, levels = str_sort(unique(mfi_causal$ST),decreasing=T))) %>%
    spread(grp,rho)%>% 
    ggplot() + 
    geom_boxplot(aes(ST,surr)) +
    geom_point(aes(x=ST,y=raw, shape=sig),size=2, color='red') +
    scale_shape_manual(values = c(1, 16)) +
    labs(x='',y='') +
    ggtitle("T") + 
    scale_y_continuous(limits = c(-0.15, 0.1)) +
    coord_flip() + mytheme2 
  
  
  p3_fluP <- mfi_causal %>% filter(tp_value==lag) %>%
    filter(dis=="fluP", plt=="o3")%>% 
    mutate(sig=as.factor(ifelse(p<=0.05, 'sig', 'non_sig'))) %>%
    mutate(ST=factor(ST, levels = str_sort(unique(mfi_causal$ST),decreasing=T))) %>%
    spread(grp,rho)%>% 
    ggplot() + 
    geom_boxplot(aes(ST,surr)) +
    geom_point(aes(x=ST,y=raw, shape=sig),size=2, color='red') +
    scale_shape_manual(values = c(1, 16)) +
    labs(x='',y='') +
    ggtitle(expression(O[3]))  + 
    scale_y_continuous(limits = c(-0.15, 0.1)) +
    coord_flip() + mytheme2 
  
  p_causal_mfi <- ggarrange(p1_fluP, p2_fluP, p3_fluP, nrow = 1, ncol = 3, 
                            widths = c(1.5, 1, 1), align = "h") %>% 
    annotate_figure(bottom  = text_grob(expression(paste("Forecast Improvement (", Delta, rho, ")")), 
                                        size = 18, face="bold"))
  
  ggsave(p_causal_mfi, file=paste0(OutputPath, 
                                   "fluP_causal_mfi_", 
                                   as.character(abs(lag)),
                                   ".tiff"), 
         width=14, height=14, units = "in", scale = 0.8)
  
  print(p_causal_mfi)
  
}

mfi_causal %>% filter(dis=="fluP", plt=="o3") %>% ungroup() %>% summarise(rho=quantile(rho, c(0,1)), prob =c(0,1))

fn_plt_mfi(0)
fn_plt_mfi(-1)
fn_plt_mfi(-2)

###### to check whether delta_rho>0 in total -- Wilcox ###### 

load('mfi_out_2021flu_o3_H.rda')

one_sample_mean=function(x){
  wilcox.test(x, mu = 0, alternative = "greater")
}

df_A2=df_123_fluP %>% group_by(tp_value,plt) %>%
  summarise(mean=mean(delta_rho),median=median(delta_rho),
            p_gt_0=one_sample_mean(delta_rho)$p.value)

df_A2 %>% arrange(plt, tp_value)  # 30 raw delta_rho--Wilcox--for single+combined predictors


###### to check whether delta_rho is different from each other in total ###### 

df_0A=df_123_fluP %>% 
  select(state,tp_value,plt,delta_rho) %>%
  spread(plt,delta_rho)

df_0B=df_0A %>%
  select(-state,-tp_value)

fn_p_less=function(data,tp_values,small_V,big_V){
  temp.data=data %>% filter(tp_value==tp_values)
  p_less=wilcox.test(temp.data[,small_V], temp.data[,big_V], paired = TRUE,
                     alternative = "less")
  data.frame(tp_values,small_V,big_V,p_less$p.value)
}

plist_p=list(data=list(df_0A),
             tp_values=-2:0,
             small_V=names(df_0B),
             big_V =names(df_0B)) %>%
  cross_df() %>%
  filter(small_V != big_V, big_V %in% c('O3','oa','ot'))

delta_p_out <- plist_p %>% pmap_df(fn_p_less)

delta_p_out

##### 2) Plotting MFI test #####

mytheme <-  theme_bw() +
  theme(axis.title = element_text(size=18),
        axis.text.x=element_text(angle = 90, size = 16, colour = "black"),
        axis.text.y=element_text(size = 16, colour = "black"),
        axis.ticks.length=unit(0.1,'cm'),
        plot.title = element_text(size=18, hjust=0.5)
  )

fn_plt_mfi_vs <- function(lag){
  
  p_mfi_vs <- df_123_fluP %>% 
    filter(tp_value==lag)%>% 
    mutate(plt=factor(plt, levels = c("AH", "T", "O3", "oa",
                                      "ot", "at", "o_t_a"))) %>%
    filter(plt!="at", plt!="o_t_a") %>%
    ggplot(aes(x=factor(plt), y=delta_rho)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = 2, color = "red") +
    labs(x =  "",
         y = expression(paste("Forecast Improvement (", Delta, rho, ")"))) +
    scale_x_discrete(labels = c(expression("AH"),
                                expression("T"),
                                expression(O[3]),
                                expression(AH+~O[3]),
                                expression(T+~O[3]))) +
    scale_y_continuous(limits=c(-0.015,0.025), breaks=seq(-0.015,0.025, by= 0.005)) +
    mytheme
  
  ggsave(p_mfi_vs, file=paste0(OutputPath, 
                               "fluP_mfi_vs_", 
                               as.character(abs(lag)),
                               ".tiff"), 
         width=7, height=7, units = "in", scale = 0.8)
  
  print(p_mfi_vs)
}

p_mfi_vs_0 <- fn_plt_mfi_vs(0)
p_mfi_vs_1 <- fn_plt_mfi_vs(-1)
p_mfi_vs_2 <- fn_plt_mfi_vs(-2)

p_mfi_vs_0_merge <- p_mfi_vs_0 + ggtitle("Lag 0")
p_mfi_vs_1_merge <- p_mfi_vs_1 + ggtitle("Lag 1") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())
p_mfi_vs_2_merge <- p_mfi_vs_2 + ggtitle("Lag 2") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())

p_mfi_vs <- ggarrange(p_mfi_vs_0_merge, p_mfi_vs_1_merge, p_mfi_vs_2_merge,
                      nrow = 1, ncol = 3,align = "h",
                      widths = c(1.4, 1, 1))
p_mfi_vs 

# ggsave(p_mfi_vs,filename = paste0(OutputPath,'p_mfi_vs.tiff'), 
#        width=9, height=7, units = "in", scale = 0.8)
# 

##### 3) Plotting O3 MFI Map #####

unique(df_123_fluP$plt)

fn_mfi_map <- function(plt, lag){
  
  mfi_fluP <- df_123_fluP %>% 
    filter(plt==plt, tp_value==lag) %>% 
    mutate(state=toupper(state)) %>%
    select(state, delta_rho)
  
  p_mfi_map <-  plot_usmap(data = mfi_fluP, values = 'delta_rho', 
                           color = "black", labels = T) + 
    scale_fill_continuous(name = expression(bold(paste(Delta,rho))),
                          low="#74add1", high="#d73027",
                          guide="colorbar", na.value="white",
    ) + 
    theme(legend.position = "right",
          legend.title=element_text(size=18),
          legend.text = element_text(size=16, color = "black"),
          legend.key.size = unit(1, 'cm'))
  
  ggsave(p_mfi_map,filename = paste0(OutputPath,'p_', as.character(plt),'_mfi_map_',
                                     as.character(abs(lag)),'.tiff'),
         width=8, height=6, units = "in", scale = 1.5)
  
  print(p_mfi_map)
  
}

fn_mfi_map("O3", 0)
fn_mfi_map("O3", -1)
fn_mfi_map("O3", -2)

df_123_fluP %>% filter(plt=="O3", tp_value==0) %>% 
  summarise(delta_rho=quantile(delta_rho, seq(0,1, 0.1), na.rm=T), prob =seq(0,1, 0.1))

############################## 3. S-map: effect strength ############################

##### 1) Prepare data for plotting #####

load('C_out_2021flu_o3_H.rda')

names(C_out)

SEeffect = C_out %>% 
  dplyr::rename(state='ST') %>%
  select(date,state,tp_value,dis,plt,effect) %>%
  mutate(plt=factor(plt,
                    levels=c("ah", "temp", "o3"),
                    labels=c("AH", "T", "O[3]"))) %>%
  mutate(tp_value=factor(tp_value,
                         levels=c(0, -1, -2),
                         labels=c('Lag 0','Lag 1','Lag 2'))) %>%
  left_join(df_T,by=c('date', 'state')) 

names(SEeffect)

SEeffect %>% 
  group_by(tp_value,plt,dis) %>%
  dplyr::summarise(n = n(),
                   qs = quantile(effect, 
                                 c(0, 0.01, 0.05, 0.5, 0.95, 0.99, 1),
                                 na.rm=T), 
                   prob = c(0, 0.01, 0.05, 0.5, 0.95, 0.99, 1)) %>%
  print(n= 100)

SEeffect_ex = SEeffect %>% group_by(tp_value,plt,dis) %>% 
  filter(effect < quantile(effect, probs=.99, na.rm = T), 
         effect > quantile(effect, probs=.01, na.rm = T),)

SEeffect_ex %>% group_by(tp_value,plt,dis) %>%
  dplyr::summarise(n = n(),
                   qs = quantile(effect, 
                                 c(0, 0.01, 0.05, 0.5, 0.95, 0.99, 1)), 
                   prob = c(0, 0.01, 0.05, 0.5, 0.95, 0.99, 1)) %>% 
  print(n= 100)

##### 2) Plotting o3 effect ABCD #####

mytheme3 <-  theme_bw() +
  theme(axis.title = element_text(size=18),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(vjust= 2.5),
        axis.text.x=element_text(size = 16, colour = "black"),
        axis.text.y=element_text(size = 16, colour = "black"),
        axis.ticks.length=unit(0.1,'cm'),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(size = rel(1.2), face = "bold"),
        plot.title = element_text(size=18, hjust=0.5),
        plot.margin = unit(c(0,0.1,0.5,0.5), "lines"))

###############

fn_plt_smap_box <- function(lag){
  
  SEeffect_o3 <- SEeffect_ex %>%
    filter(dis=="fluP", plt=="O[3]") %>%
    filter(tp_value==lag)
  
  p_o3effect_t <- ggplot(SEeffect_o3, 
                         aes(date, effect)) +
    geom_point(color = "grey", alpha = 0.2, size = 1.5) +
    geom_hline(yintercept = 0, linetype = 2, color = "red") +
    labs(x = " ", 
         y =  expression(atop("Effect of "* O[3]* " on Influenza", "Influenza/ "* O[3]))) +
    scale_y_continuous(limits=c(-0.40, 0.10), 
                       breaks = seq(-0.40, 0.10, by=0.10),
                       labels = scales::number_format(accuracy = 0.01)) +
    scale_x_date(breaks = seq(as.Date("2011-01-01"), 
                              as.Date("2015-01-01"), by="1 year"), 
                 date_labels="%Y") + mytheme3 +
    theme(axis.title.x =element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x = element_blank())
  
  p_fluP_t <- ggplot(SEeffect_o3, aes(date, fluP)) +
    geom_point(color = "grey", alpha = 0.2, size = 1.5) +
    labs(x = "Date", 
         y = expression(atop("Influenza Activity", ""))) +
    scale_y_continuous(limits=c(0, 0.09), breaks = seq(0, 0.09, by=0.02),
                       labels = scales::number_format(accuracy = 0.01))+
    scale_x_date(breaks = seq(as.Date("2011-01-01"), 
                              as.Date("2015-01-01"), by="1 year"), 
                 date_labels="%Y") + mytheme3 
  
  SE_o3_fluP_month = SEeffect_o3 %>%
    mutate(month=month(ymd(date), label = TRUE, abbr = TRUE)) 
  
  p_o3effect_box <- SE_o3_fluP_month %>%
    ggplot(aes(month, effect)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = 2, color = "red") +
    scale_y_continuous(limits=c(-0.40, 0.10), 
                       breaks = seq(-0.40, 0.10, by=0.10),
                       labels = scales::number_format(accuracy = 0.01))+
    scale_x_discrete(limits=c('Sep','Oct','Nov',
                              'Dec','Jan','Feb',
                              'Mar','Apr','May',
                              'Jun','Jul','Aug')) +
    labs(x = " ", 
         y = " ") + mytheme3 +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())
  
  p_fluP_box <- SE_o3_fluP_month %>%
    ggplot(aes(month, fluP)) +
    geom_boxplot() +
    theme_bw() +
    scale_y_continuous(limits=c(0.00, 0.09), breaks = seq(0.00, 0.09, by=0.02),
                       labels = scales::number_format(accuracy = 0.01))+
    scale_x_discrete(limits=c('Sep','Oct','Nov',
                              'Dec','Jan','Feb',
                              'Mar','Apr','May',
                              'Jun','Jul','Aug')) +
    labs(x = "Month", 
         y = " ") + mytheme3 +
    theme(axis.text.y=element_blank(),
          axis.ticks.y = element_blank())
  
  o3effects_fluP_box <- ggarrange(p_o3effect_t, 
                                  p_o3effect_box,
                                  p_fluP_t, 
                                  p_fluP_box,
                                  labels = c("A", "B", "C", "D"),hjust=-0.8,
                                  font.label = list(size = 20, color = "black", face = "bold"),
                                  ncol = 2, nrow = 2, align = "hv")
  
  ggsave(o3effects_fluP_box, 
         filename = paste0(OutputPath,
                           'o3effects_fluP_box',
                           substr(lag, 5,5),'.tiff'),
         width=14, height=12, units = "in", scale = 0.9)
  
  print(o3effects_fluP_box)
}

fn_plt_smap_box("Lag 0")
fn_plt_smap_box('Lag 1')
fn_plt_smap_box('Lag 2')

####### HI and AK separately #######

o3_effect_hi <- ggplot(SEeffect_ex %>%
                         filter(dis=="fluP", plt=="O[3]") %>% 
                         filter(state == c('hi')), 
                       aes(date, effect)) +
  facet_grid(tp_value ~ ., scales = 'free')+
  geom_point(color = "grey", alpha = 0.8, size = 1.5) +
  geom_hline(yintercept = 0, linetype = 2, color = "red") +
  labs(x = " ", 
       y =  expression(atop("Effect of "* O[3]* " on Influenza", "Influenza/ "* O[3]))) +
  ggtitle('Hawaii') +
  scale_y_continuous(limits=c(0, 0.10), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_x_date(breaks = seq(as.Date("2011-01-01"), 
                            as.Date("2015-01-01"), by="1 year"), 
               date_labels="%Y") + mytheme3

o3_effect_hi

# ggsave(o3_effect_hi,filename = paste0(OutputPath,"/excludeHiAk/",'o3_effect_hi.tiff'),
#        width=8, height=8, units = "in", scale = 0.8)


o3_effect_ak <-ggplot(SEeffect_ex %>%
                        filter(dis=="fluP", plt=="O[3]") %>%
                        filter(state == c('ak')), 
                      aes(date, effect)) +
  facet_grid(tp_value ~ .)+
  geom_point(color = "grey", alpha = 0.8, size = 1.5) +
  geom_hline(yintercept = 0, linetype = 2, color = "red") +
  labs(x = " ", 
       y =  expression(atop("Effect of "* O[3]* " on Influenza", "Influenza/ "* O[3]))) +
  ggtitle('Alaska') +
  scale_y_continuous(limits=c(-0.05, 0.05), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_x_date(breaks = seq(as.Date("2011-01-01"), 
                            as.Date("2015-01-01"), by="1 year"), 
               date_labels="%Y") + mytheme3
o3_effect_ak

# ggsave(o3_effect_ak,filename = paste0(OutputPath,"/excludeHiAk/",'o3_effect_ak.tiff'),
#        width=8, height=8, units = "in", scale = 0.8)

##### 3) correlation coefficient calculation #####

cor.mat <- matrix(nrow=3, ncol=2, dimnames=list(c("lag0", "lag1","lag2"),
                                                c("overall","month_mean")))

for (i in 1:3) {
  SEeffect_o3 <- SEeffect_ex %>%
    filter(dis=="fluP", plt=="O[3]") %>%
    filter(tp_value==rev(unique(SEeffect_ex$tp_value))[i])
  
  SE_o3_fluP_month = SEeffect_o3 %>%
    mutate(month=month(ymd(date), label = TRUE, abbr = TRUE))
  
  o3effect_fluP_month_mean = SE_o3_fluP_month %>%  
    select(effect, fluP, month) %>% 
    group_by(month) %>% 
    dplyr::summarise(effect=mean(effect,na.rm=TRUE),
                     fluP=mean(fluP,na.rm=TRUE))
  
  cor.mat[i,1] = cor.test(SEeffect_o3$effect, SEeffect_o3$fluP)$estimate
  
  cor.mat[i,2] = cor.test(o3effect_fluP_month_mean$effect, 
                          o3effect_fluP_month_mean$fluP)$estimate
}
print(cor.mat)

p.mat <- matrix(nrow=3, ncol=2, dimnames=list(c("lag0", "lag1","lag2"),
                                              c("overall","month_mean")))

for (i in 1:3) {
  SEeffect_o3 <- SEeffect_ex %>%
    filter(dis=="fluP", plt=="O[3]") %>%
    filter(tp_value==rev(unique(SEeffect_ex$tp_value))[i])
  
  SE_o3_fluP_month = SEeffect_o3 %>%
    mutate(month=month(ymd(date), label = TRUE, abbr = TRUE))
  
  o3effect_fluP_month_mean = SE_o3_fluP_month %>%  
    select(effect, fluP, month) %>% 
    group_by(month) %>% 
    dplyr::summarise(effect=mean(effect,na.rm=TRUE),
                     fluP=mean(fluP,na.rm=TRUE))
  
  p.mat[i,1] = cor.test(SEeffect_o3$effect, SEeffect_o3$fluP)$p.value
  
  p.mat[i,2] = cor.test(o3effect_fluP_month_mean$effect, 
                        o3effect_fluP_month_mean$fluP)$p.value
}
print(p.mat)

##### 4) Plotting AH T effect by date #####

p_effect_t <- SEeffect_ex %>%
  filter(dis=="fluP")  %>%
  filter(plt!="O[3]")  %>%
  ggplot(aes(date, effect)) +
  geom_point(color = "grey", alpha = 0.2, size = 1.5) +
  facet_grid(tp_value ~ plt, scales = 'free') +
  geom_hline(yintercept = 0, linetype = 2, color = "red") +
  theme_bw() +
  labs(x = "Date", y = "Effect on Influenza") +
  scale_x_date(breaks = seq(as.Date("2011-01-01"),
                            as.Date("2015-01-01"), by="1 year"),
               date_labels="%Y") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+ 
  theme_bw() +
  theme(axis.title = element_text(size=18),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(vjust= 2.5),
        axis.text.x=element_text(size = 12, colour = "black"),
        axis.text.y=element_text(size = 12, colour = "black"),
        axis.ticks.length=unit(0.1,'cm'),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size=18, hjust=0.5),
        plot.margin = unit(c(0,0.1,0.5,0.5), "lines"))+
  theme(strip.text = element_text(size = rel(1.3), face = "bold"))   

p_effect_t

# ggsave(p_effect_t,filename = paste0(OutputPath, 'p_effect_t.tiff'),
#        width=10, height=8, units = "in", scale = 0.8)

fn_plt_effect_t <- function(lag){
  
  p_effect_t <- SEeffect_ex %>%
    filter(dis=="fluP")  %>%
    filter(plt!="O[3]")  %>%
    filter(tp_value==lag) %>%
    ggplot(aes(date, effect)) +
    geom_point(color = "grey", alpha = 0.4, size = 1.5) +
    facet_grid(plt ~ ., scales = 'free', labeller = label_parsed) +
    geom_hline(yintercept = 0, linetype = 2, color = "red") +
    theme_bw() +
    labs(x = "Date", y = "Effect on Influenza") +
    scale_x_date(breaks = seq(as.Date("2011-01-01"),
                              as.Date("2015-01-01"), by="1 year"),
                 date_labels="%Y") +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+ 
    mytheme3 
  
  ggsave(p_effect_t, 
         filename = paste0(OutputPath,
                           'p_effect_t_',
                           substr(lag, 5,5),'.tiff'),
         width=8, height=7, units = "in", scale = 0.8)
  
  print(p_effect_t)
}

# fn_plt_effect_t("Lag 0")
# fn_plt_effect_t("Lag 1")
# fn_plt_effect_t("Lag 2")


##### 5) Plotting O3 effect size Map #####

fn_SE_map <- function(lag){
  O3_SE_mean = SEeffect_ex %>%
    filter(dis=="fluP", plt=="O[3]") %>%
    filter(tp_value==lag) %>%
    group_by(state) %>% 
    dplyr::summarise(o3_effect=mean(effect,na.rm=TRUE)) %>%
    mutate(state=toupper(state))
  
  range(O3_SE_mean$o3_effect)
  
  p_SE_o3map <-  plot_usmap(data = O3_SE_mean, values = 'o3_effect', 
                            color = "black", labels = T) + 
    scale_fill_continuous(name = "Effect",
                          limits = c(-0.10,0.05), breaks = c(-0.10, -0.05, 0, 0.05),
                          low="#3288bd", high="#d53e4f",
                          guide="colorbar", na.value="white") + 
    theme(legend.position = "right",
          legend.title=element_text(size=18),
          legend.text = element_text(size=16, color = "black"),
          legend.key.size = unit(1, 'cm'))
  
  ggsave(p_SE_o3map,filename = paste0(OutputPath,'p_SE_o3map_',
                                      substr(lag, 5,5),'.tiff'),
         width=8, height=6, units = "in", scale = 1.5)
  
  print(p_SE_o3map)
}

fn_SE_map("Lag 0")
fn_SE_map("Lag 1")
fn_SE_map("Lag 2")

