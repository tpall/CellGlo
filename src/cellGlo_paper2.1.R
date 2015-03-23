library('ProjectTemplate')
rm(list=ls())
load.project()

# Two objectives:
#1. lets try to normalise 3MUT-treated data against control treatment
#2. create plate mask of ege and region effects

tec %>% summary

# Remove edge effect stuff
# First in 140610 experiment
tec$value[tec$date=="140610"] <- tec %>%
  filter(!treatment=="media"&date == "140610") %>%
  with(tapply(value,list(row,col),function(x) x)) %>%
  medpolish(na.rm = TRUE) %>%
  with((overall + outer(rep(1,8),col, "+") + residuals)) %>% t %>% c 

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

tec %>% 
  filter(treatment%in%c("rhIgG-Fc","3MUT-Fc")) %>%
  group_by(date,doses,GF) %>%
  mutate(value=value/median(value[treatment=="rhIgG-Fc"],na.rm = TRUE)) %>%
  dlply(.(GF),Plotfun) %>%
  do.call(grid.arrange,.)
