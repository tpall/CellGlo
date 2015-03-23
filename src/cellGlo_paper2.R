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
# ivi <- snc %>% filter(grepl("IVIS",exp.id))

# add variables content for normalisation and treat2
tec %<>%
  mutate(content = ifelse(treatment=="UT","pos","sample"),
         content = ifelse(doses_GF==0,"neg",content),
         content = ifelse(treatment=="media","blank",content),
         treat2 = paste(doses_GF,GF,treatment))

# Before normalisation
tec %>% 
  group_by(date,plate) %>% 
#   mutate(value=scale(value)) %>%
  ggplot(aes(x=factor(row),y=value,shape=factor(GF),color=treatment))+
  geom_point(position = position_dodge(0.5)) +
  facet_wrap(~date,scale="free_y") +
  scale_color_colorblind()

# Run medpolish to remove row effect and insert back corrected data ----
tec$value[tec$date=="140610"] <- tec %>%
  filter(!treatment=="media"&date == "140610") %>%
  with(tapply(value,list(row,col),function(x) x)) %>%
  medpolish(na.rm = TRUE) %>%
  with((overall + outer(rep(1,8),col, "+") + residuals)) %>% t %>% c 

# # after norm
# tec %>%
#   filter(!treatment=="media") %>%
#   group_by(date,plate,GF) %>% 
#   mutate(value=scale(value)) %>% 
#   dlply(.(GF)) %>%
#   lapply({.%>%{
#     ggplot(.,aes(x=log10(doses),y=value,color=treat2)) +
#       stat_summary(fun.data = mean_se, geom = "pointrange") +
#       stat_summary(fun.y = mean, geom = "line")+
#       scale_color_colorblind()}})
# 
# tec %>%
#   filter(!treatment=="media") %>%
#   group_by(date,plate,GF) %>% 
#   mutate(value=scale(value)) %>% 
#   ggplot(aes(x=log10(doses),y=value,color=treat2)) +
#   stat_summary(fun.data = mean_se, geom = "pointrange") +
#   stat_summary(fun.y = mean, geom = "line") +
#   facet_grid(~GF)

# Adjust edge effect in HGF data 
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
  lapply({.%>%with(overall+outer(row,col,"+")+residuals)%>%t%>%c})

# Insert back corrected data ----
tec$value[tec$date=="140705"&tec$GF=="HGF"] <- med[[1]] 
tec$value[tec$date=="140903"&tec$GF=="HGF"] <- med[[2]] 
tec$value[tec$date=="140917"&tec$GF=="HGF"] <- med[[3]] 

med <- ivi %>%
  filter(!treatment=="media"&GF=="HGF") %>%
  group_by(date,plate) %>%
  dlply(.(date)) %>%
  lapply({.%>%with(tapply(value,list(row,col),function(x) x)) %>% medpolish(na.rm = TRUE)}) %>% 
  lapply({.%>%{
    out <- Predfun(.)
    Medfun(.,out)
  }}) %>%
  lapply({.%>%with(overall+outer(row,col,"+")+residuals)%>%t%>%c})

# Insert back corrected data ----
ivi$value[ivi$date=="140705"&ivi$GF=="HGF"] <- med[[1]] 
ivi$value[ivi$date=="140903"&ivi$GF=="HGF"] <- med[[2]] 
ivi$value[ivi$date=="140917"&ivi$GF=="HGF"] <- med[[3]] 

# Plot summary results
# Original plot
tec %>%
  filter(!treatment=="media") %>%
  group_by(date,plate,GF) %>% 
  mutate(value=scale(value)) %>% 
  ggplot(aes(x=log10(doses),y=value,color=treatment,shape=factor(doses_GF))) +
  stat_summary(fun.data = mean_se, geom = "pointrange") +
  stat_summary(fun.y = mean, geom = "line")+
  facet_wrap(~GF,scale="free_y")+
  scale_color_colorblind() 

Plotfun <- .%>% {
  ggplot(.,aes(x=log10(doses),y=value,color=treatment)) +
    stat_summary(fun.data = mean_se, geom = "pointrange") +
        stat_summary(fun.y = mean, geom = "line",size=1)+
    scale_color_colorblind()}

ntec <- tec %>%
  filter(!treatment=="media"&treatment%in%c("3MUT-Fc","rhIgG-Fc")) %>%
  group_by(date,plate) %>% 
  mutate(value=scale(value)) %>%
  group_by(date,plate,GF,doses,doses_GF,treatment,treat2) %>% 
  summarise(value=mean(value))

ntec %>% 
  dlply(.(GF),Plotfun) %>%
  do.call(grid.arrange,.)

ntec2 <- tec %>%
  filter(!treatment=="media") %>%
  group_by(date,plate,GF) %>% 
  mutate(value=(value-mean(value[content=="neg"]))/(mean(value[content=="pos"])-mean(value[content=="neg"]))) %>%
  group_by(date,plate,GF,doses,doses_GF,treatment,treat2) %>% 
  summarise(value=mean(value))

ntec2 %>% 
  filter(treatment%in%c("3MUT-Fc","rhIgG-Fc")) %>%
  dlply(.(GF),Plotfun) %>%
  do.call(grid.arrange,.)