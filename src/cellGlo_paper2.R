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

# median normalisation
p <- tec %>%
  mutate(content = ifelse(treatment=="UT","pos","sample"),
         content = ifelse(doses_GF==0,"neg",content),
         content = ifelse(treatment=="media","blank",content)) %>% 
  group_by(date,plate) %>% 
  mutate(value=value/median(value[content=="sample"])) %>% # apply median filter
  mutate(value=(value-mean(value[content=="sample"]))/sd(value[content=="sample"])) %>% # standardize
  mutate(treat2 = paste(doses_GF,GF,treatment)) %>%
  group_by(date,plate,GF,doses,treatment,treat2) %>%
  summarise(value=sqrt(sum(value^2)/length(value))) %>% # RMS mean
  filter(!treatment=="media") %>%
  dlply(.(date)) %>%
  lapply({.%>%{
    ggplot(.,aes(x=log10(doses),y=value,color=treat2)) +
      stat_summary(fun.data = mean_sdl, mult=1, geom = "pointrange") +
      stat_summary(fun.y = mean, geom = "line") +
      facet_grid(~GF,scale = "free_y")}})  
  
  summary
  mutate(content=as.factor(content)) %>% 
  
  
  p <- tec %>% 
  mutate(content = ifelse(treatment=="UT","pos","sample"),
         content = ifelse(doses_GF==0,"neg",content),
         content = ifelse(treatment=="media","blank",content)) %>% 
  filter(!treatment=="media") %>%
  group_by(date,plate) %>% 
  mutate(value=value/median(value[content=="sample"])) %>% # apply median filter
  mutate(value=(value-mean(value[content=="sample"]))/sd(value[content=="sample"])) %>% # standardize
  mutate(treat2 = paste(doses_GF,GF,treatment)) %>%
  group_by(date,plate,GF,doses,treatment,treat2,content) %>%
  summarise(value=sqrt(sum(value^2)/length(value))) %>% # RMS mean
  filter(content=="sample") %>%
  dlply(.(GF)) %>%
  lapply({.%>%{
    ggplot(.,aes(x=log10(doses),y=value,color=treat2)) +
      stat_summary(fun.data = mean_se, geom = "pointrange") +
      stat_summary(fun.y = mean, geom = "line") +
      scale_color_colorblind()}})

grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],ncol=2)
