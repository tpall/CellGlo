library('ProjectTemplate')
rm(list=ls())
load.project()

# query list object col numbers to find "Sync" experiments;snc=12,async=10,sirna=9
snc <- datalist %>% 
  lapply({.%>%names%>%length}) %>% 
  unlist %>% 
  equals(12) %>% 
  datalist[.] %>%
  bind_rows

tec <- snc %>%
  filter(!grepl("IVIS",exp.id))

tec %>% 
  group_by(date,plate) %>%
  mutate(value=scale(value,center = FALSE)%>%c)%>%
  filter(treatment=="media") %>%
  summarise(Mean=mean(value),
            SD=sd(value),
            Median=median(value),
            Max=max(value))

tec %>%
  group_by(date,plate) %>%
  mutate(value = value - mean(value[treatment == "media"], na.rm = TRUE)) %>%
  mutate(value = scale(value,center = FALSE)) %>%
  filter(!treatment=="media") %>%
  qplot(x=factor(date),y=value,geom="boxplot",data=.,fill=factor(plate))

# log transform values
df$doses <- (df$doses + 0.1)/1e9
