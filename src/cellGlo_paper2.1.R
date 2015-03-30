library('ProjectTemplate')
rm(list=ls())
load.project()

# Local munge
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
         content = ifelse(treatment=="FUM","xFUM",content),
         treat2 = paste(doses_GF,GF,treatment))

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
    stat_summary(fun.y = mean, geom = "bar") +
    scale_x_discrete(breaks=c("neg","pos","xFUM"),
                     labels=c("Unind.", GF, paste(GF,"+\nFUM"))) +
    xlab(NULL) +
    ylab("Relative cell number") +
    theme_classic() +
    theme(axis.text.x=element_text(angle = 90)) +
    scale_y_continuous(expand = c(0, 0)) +
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
#   GF <- select(.,GF) %>% "["(1,1) %>% as.character %>% ifelse(.=="bFGF","FGF2",.) %>% as.character
  p <- ggplot(.,aes(x=log10(doses),y=value,shape=treatment)) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width=0.075) +
    stat_summary(fun.y = mean, geom = "point",size=3) +
    stat_summary(fun.y = mean, geom = "line") +
    ylab(NULL) +
    xlab(bquote(list(Conc.,log[10](M)))) +
    scale_y_continuous(expand = c(0, 0)) +
    expand_limits(y = 0)+
    theme_classic() +
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
labs <- c("VEGF","GDF-2","HGF","FGF2") %>% lapply(textGrob)

fig <- do.call(arrangeGrob, c(limp, list(ncol=4,widths=c(2,3))))
fig

# lets fit some models
# tec %>%
#   filter(!treatment=="media") %>% # empty wells
#   filter(doses<1.50001e-05) %>%
#   group_by(date,plate,GF) %>% 
#   mutate(value=(value-min(value))/(range(value)%>%diff)) %>%
#   ggplot(aes(x=log10(doses),y=value,color=treatment,shape=factor(doses_GF))) +
#   stat_summary(fun.data = mean_se, geom = "pointrange") +
#   stat_summary(fun.y = mean, geom = "line")+
#   facet_wrap(~GF)+
#   scale_color_colorblind()