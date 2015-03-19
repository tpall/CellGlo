library('ProjectTemplate')
rm(list=ls())
load.project()

library(ggthemes)

# query list object col numbers to find "Sync" experiments;snc=12,async=10,sirna=9 ----
snc <- datalist %>% 
  lapply({.%>%names%>%length}) %>% 
  unlist %>% 
  equals(12) %>% 
  datalist[.] %>%
  bind_rows

snc %<>% mutate(doses = doses%>%as.numeric,
                doses_GF = doses_GF%>%as.numeric,
                GF = GF%>%as.factor,
                date = date%>%as.factor)

# log transform treatment doses ----
snc %<>% mutate(doses = (doses+0.1)/1e9)

# lets use only Tecan measurements ----
tec <- snc %>% filter(!grepl("IVIS",exp.id))
ivi <- snc %>% filter(grepl("IVIS",exp.id))

# add variables content for normalisation and treat2
tec %<>%
  mutate(content = ifelse(treatment=="UT","pos","sample"),
         content = ifelse(doses_GF==0,"neg",content),
         content = ifelse(treatment=="media","blank",content),
         treat2 = paste(doses_GF,GF,treatment))

# befor norm
tec %>% 
  group_by(date,plate) %>% 
  mutate(value=scale(value)) %>%
  ggplot(aes(x=factor(row),y=value,shape=factor(GF),color=treatment))+
  geom_point(position = position_dodge(0.5)) +
  facet_wrap(~date,scale="free_y") +
  scale_color_colorblind()

# We have strong row effect in 140610 experiment
# Lets try to get rid of row effect by median polish
med <- tec %>%
  filter(!treatment=="media"&date == "140610") %>%
  with(tapply(value,list(row,col),function(x) x)) %>%
  medpolish(na.rm = TRUE)

# Plot results of median polish
(med$overall + outer(rep(1,8),med$col, "+") + med$residuals) %>% 
  melt %>% 
  set_colnames(c("row","col","value")) %>%
  ggplot(aes(x=factor(row),y=value)) +
  geom_point() 

# Insert back corrected data ----
tec$value[tec$date=="140610"] <- (med$overall+outer(rep(1,8),med$col,"+")+med$residuals)%>%t%>%c 

# after norm
tec %>%
  filter(!treatment=="media") %>%
  group_by(date,plate,GF) %>% 
  mutate(value=scale(value)) %>% 
  dlply(.(GF)) %>%
  lapply({.%>%{
    ggplot(.,aes(x=log10(doses),y=value,color=treat2)) +
      stat_summary(fun.data = mean_se, geom = "pointrange") +
      stat_summary(fun.y = mean, geom = "line")+
      scale_color_colorblind()}})
  
# HGF data 
tec %>%
  filter(!treatment=="media"&GF=="HGF") %>%
  group_by(date,plate,GF) %>% 
  mutate(value=scale(value)) %>% 
  ggplot(.,aes(x=log10(doses),y=value,color=treat2)) +
  stat_summary(fun.data = mean_se, geom = "pointrange") +
  stat_summary(fun.y = mean, geom = "line")+
#   facet_wrap(~date) +
  scale_color_colorblind()

tec %>%
  filter(!treatment=="media"&GF=="HGF") %>%
  group_by(date,plate,GF) %>% 
  mutate(value=scale(value)) %>% 
  ggplot(.,aes(x=row,y=value,color=treat2)) +
  stat_summary(fun.data = mean_se, geom = "pointrange") +
  stat_summary(fun.y = mean, geom = "line")+
  facet_wrap(~date) +
  scale_color_colorblind()


# 
# HGF data 
Predfun <- function(x) x$col[2:3]%>%
  data.frame(col=c(2,3),eff=.)%>%
  lm(eff~col,data=.)%>%
  predict(.,data.frame(col=1))

Medfun <- function(x,newfirst) {
  x$col[1]<-newfirst
  x}

med <- tec %>%
  filter(!treatment=="media"&GF=="HGF") %>%
  group_by(date,plate) %>%
  dlply(.(date)) %>%
  lapply({.%>%with(tapply(value,list(row,col),function(x) x)) %>% medpolish(na.rm = TRUE)}) %>% 
  lapply({.%>%{
    out <- Predfun(.)
    Medfun(.,out)
  }}) %>%
  lapply(function(x) (x$overall+outer(x$row,x$col,"+")+x$residuals)%>%t%>%c)

# Insert back corrected data ----
tec$value[tec$date=="140705"&tec$GF=="HGF"] <- med[[1]] 
tec$value[tec$date=="140903"&tec$GF=="HGF"] <- med[[2]] 
tec$value[tec$date=="140917"&tec$GF=="HGF"] <- med[[3]] 

# p <- tec %>% 
#   filter(!treatment=="media") %>%
#   group_by(date,plate) %>% 
#   mutate(value=value/median(value[content=="sample"])) %>% # apply median filter
#   mutate(value=(value-mean(value[content=="sample"]))/sd(value[content=="sample"])) %>% # standardize
#   group_by(date,plate,GF,doses,treatment,treat2,content) %>%
#   summarise(value=sqrt(sum(value^2)/length(value))) %>% # RMS mean
#   filter(content=="sample") %>%
#   dlply(.(GF)) %>%
#   lapply({.%>%{
#     ggplot(.,aes(x=log10(doses),y=value,color=treat2)) +
#       stat_summary(fun.data = mean_se, geom = "pointrange") +
#       stat_summary(fun.y = mean, geom = "line") +
#       scale_color_colorblind()}})
# 
# grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],ncol=2)
