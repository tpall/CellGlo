library('ProjectTemplate')
rm(list=ls())
load.project()

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

# about background ----
tec %>% 
  group_by(date,plate) %>%
  mutate(value=scale(value,center = FALSE)%>%c)%>%
  filter(treatment=="media") %>%
  summarise(Mean=mean(value),
            SD=sd(value),
            Median=median(value),
            Max=max(value))
tec %>%
  mutate(back = ifelse(treatment=="media","bckg","cells")) %>%
  ggplot(aes(x=factor(date),y=value,fill=factor(plate),color=back)) +
  geom_boxplot() +
  stat_summary(fun.y=max,geom="point",shape=3,size=5,position = position_dodge(0.9))

# lets throw out media values from 140917 experiment plate 2, because these values seem too high ----
tec %<>% filter(!(treatment=="media"&plate=="2"&date=="140917")) 

# subtract mean background
tec %<>%
  group_by(date) %>%
  mutate(value = value - mean(value[treatment == "media"], na.rm = TRUE)) %>% 
  filter(!treatment == "media")

# lets scale values
tec %<>%
  group_by(date,plate) %>% 
  mutate(norm.value = scale(value, center = FALSE)%>%c) 

tec %>%
  ggplot(aes(x=date,y=norm.value,fill=factor(plate))) +
  geom_boxplot()

# summarise data ----
tec %>% summary

Mynorm <- function(x) (x-min(x))/(max(x)-min(x))

ivi %>%
  filter(!treatment == "media") %>%
  mutate(treat2 = paste(doses_GF,GF,treatment)) %>% # doses_GF,GF,treatment
  group_by(date,GF) %>%
  mutate(value=scale(value,center=FALSE)) %>%
  group_by(date,GF,doses,treat2) %>%
  summarise(Mean = mean(value),
            SD = sd(value)) %>%
  ggplot(aes(x=log10(doses),y=Mean,color=treat2)) +
  stat_summary(fun.data = mean_sdl, mult=1, geom = "pointrange") +
  stat_summary(fun.y = mean, geom = "line") +
  facet_wrap(~GF, scale="free_y") +
  guides(color=guide_legend(ncol=2))

library(ggthemes)
p <- ivi %>%
  filter(!treatment == "media") %>%
  mutate(treat2 = paste(doses_GF,GF,treatment)) %>% # doses_GF,GF,treatment
  group_by(date,GF) %>%
  mutate(value=scale(value)) %>%
  dlply(.(date,GF)) %>%
  lapply({.%>%{
    ggplot(.,aes(x=log10(doses),y=value,color=treat2)) +
      stat_summary(fun.data = mean_sdl, mult=1, geom = "pointrange") +
      stat_summary(fun.y = mean, geom = "line") +
      scale_color_colorblind()}})

library(gridExtra)
grid.arrange(p[[1]],p[[2]])
grid.arrange(p[[3]],p[[4]],p[[5]],p[[6]],nrow=2)
grid.arrange(p[[7]],p[[8]],p[[9]],p[[10]],nrow=2)
grid.arrange(p[[11]],p[[12]],p[[13]],p[[14]],nrow=2)

mp <- tec %>%
  filter(!treatment == "media") %>%
  dlply(.(date,plate)) %>%
  lapply({.%>%with(tapply(value,list(row,col),function(x)x))}) %>%
  lapply(medpolish,na.rm=T)

polished <- mp %>% 
  lapply({.%>%"[["(4)}) %>% 
  lapply({.%>%t%>%c%>%data.frame(polished=.)}) %>%
  bind_rows

p <- cbind(tec,polished) %>%
  filter(!treatment == "media") %>%
  mutate(treat2 = paste(doses_GF,GF,treatment)) %>% # doses_GF,GF,treatment
  dlply(.(date,GF)) %>%
  lapply({.%>%{
    ggplot(.,aes(x=log10(doses),y=polished,color=treat2)) +
      stat_summary(fun.data = mean_sdl, mult=1, geom = "pointrange") +
      stat_summary(fun.y = mean, geom = "line") +
      scale_color_colorblind()}})

library(gplots)
ivi %>%
  filter(!treatment == "media") %>%
  dlply(.(date,plate)) %>%
  lapply({.%>%with(tapply(value,list(row,col),function(x)x))}) %>%
  lapply(heatmap.2,scale = "column")