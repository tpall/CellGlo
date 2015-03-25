library('ProjectTemplate')
rm(list=ls())
load.project()

# Two objectives:
#1. lets try to normalise 3MUT-treated data against control treatment
#2. create plate mask of ege and region effects

tec %>% summary

# Remove edge effect stuff
# First in 140610 experiment
# tec$value[tec$date=="140610"] <- tec %>%
#   filter(!treatment=="media"&date == "140610") %>%
#   with(tapply(value,list(row,col),function(x) x)) %>%
#   medpolish(na.rm = TRUE) %>%
#   with((overall + outer(rep(1,8),col, "+") + residuals)) %>% t %>% c 

# HGF data are also way off 
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

library(scales)
# Stimulation plots
Plotfun <- .%>% {
  GF <- select(.,GF) %>% "["(1,1) %>% as.character %>% ifelse(.=="bFGF","FGF2",.) %>% as.character
  p <- ggplot(.,aes(x=factor(content),y=value)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width=0.25) +
  stat_summary(fun.y = mean, geom = "point",size=3) +
  scale_x_discrete(breaks=c("neg","pos","xFUM"),
                   labels=c("unind.", "ind.", "ind.+\n20 nM FUM")) +
  xlab(" ") +
  ylab("Relative cell number") +
  expand_limits(y = 0) 
# p <- p + ggtitle(paste(GF))
p}

stimp <- tec %>%
  filter(!treatment=="media") %>%
  filter(doses<1.50001e-05) %>%
  group_by(date,plate,GF) %>% 
  mutate(value=(value-min(value))/(range(value)%>%diff)) %>%
  filter(content%in%c("pos","neg","xFUM")) %>%
  mutate(value=value/median(value[content=="pos"],na.rm = TRUE)) %>%
  dlply(.(GF)) %>%
  lapply({.%>%Plotfun})

# Inhibition plots
Moreplotfun <- .%>% {
  GF <- select(.,GF) %>% "["(1,1) %>% as.character %>% ifelse(.=="bFGF","FGF2",.) %>% as.character
  p <- ggplot(.,aes(x=log10(doses),y=value,shape=treatment)) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width=0.075) +
    stat_summary(fun.y = mean, geom = "point",size=3) +
    ylab(NULL) +
    expand_limits(y = 0) +
    theme(legend.position="none") 
  # p <- p + ggtitle(paste(GF))
  p}

inhpimp <- tec %>%
  filter(!treatment=="media") %>%
  filter(doses<1.50001e-05) %>%
  group_by(date,plate,GF) %>% 
  mutate(value=(value-min(value))/(range(value)%>%diff)) %>%
  filter(treatment%in%c("rhIgG-Fc","3MUT-Fc")) %>%
  group_by(date,doses,GF) %>%
  mutate(value=value/median(value[treatment=="rhIgG-Fc"],na.rm = TRUE)) %>%
  dlply(.(GF)) %>%
  lapply({.%>%Moreplotfun})

limp <- mapply(list,stimp,inhpimp)
do.call(grid.arrange, c(limp, list(ncol=4)))

# lets fit some models
tec %>%
  filter(!treatment=="media") %>% # empty wells
  filter(doses<1.50001e-05) %>%
  group_by(date,plate,GF) %>% 
  mutate(value=(value-min(value))/(range(value)%>%diff)) %>%
  ggplot(aes(x=log10(doses),y=value,color=treatment,shape=factor(doses_GF))) +
  stat_summary(fun.data = mean_se, geom = "pointrange") +
  stat_summary(fun.y = mean, geom = "line")+
  facet_wrap(~GF)+
  scale_color_colorblind()