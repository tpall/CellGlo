tec %>% 
  ggplot(aes(x=factor(treatment),y=value,fill=factor(plate)))+
  geom_boxplot() +
  facet_wrap(~date) 

tec %>%
  filter(!(treatment=="media"&plate==2&date=="140917")) %>%
  group_by(date,plate) %>%
  mutate(value=value-(max(value[treatment=="media"]))) %>%
  filter(!treatment=="media") %>%
  ggplot(aes(x=factor(treatment),y=value,fill=factor(plate)))+
  geom_boxplot() +
  facet_wrap(~date) 

# Kvantiil normaliseerime need andmed
library(preprocessCore)
qntec <- tec %>%
  filter(!(treatment=="media"&plate==2&date=="140917")) %>%
  group_by(date) %>%
  mutate(value=value-(max(value[treatment=="media"]))) %>%
  filter(!treatment=="media") %>%
  ungroup %>%
  mutate(exp.id=paste(date,plate,sep="_"),
         well=paste(row,col,sep="")) %>% { 
           qnvals <- select(.,exp.id,well,value) %>%
             dcast(well~exp.id) %>% {
               well <- use_series(.,"well")
               cln <- select(.,-well) %>% colnames
               select(.,-well) %>% as.matrix %>% 
                 normalize.quantiles %>% data.frame %>% 
                 set_colnames(cln) %>% cbind(well,.)
             } %>%
             melt(value.name = "qnval",variable.name = "exp.id")
           left_join(.,qnvals)}

qntec %>%
  ggplot(aes(x=exp.id,y=qnval))+
  geom_boxplot()

qntec %>%
  ggplot(aes(x=log10(doses),y=qnval,color=treatment,shape=factor(doses_GF))) +
  stat_summary(fun.data = mean_se, geom = "pointrange") +
  stat_summary(fun.y = mean, geom = "line")+
  facet_wrap(~GF,scale="free_y")+
  scale_color_colorblind()

qntec %>%
  filter(!treatment=="media") %>%
  filter(doses<1.50001e-05) %>%
  group_by(date,plate,GF) %>% 
  mutate(qnval=(qnval-min(qnval))/(range(qnval)%>%diff)) %>%
  filter(treatment%in%c("rhIgG-Fc","3MUT-Fc")) %>%
  group_by(date,doses,GF) %>%
  mutate(qnval=qnval/median(qnval[treatment=="rhIgG-Fc"],na.rm = TRUE)) %>%
  ggplot(aes(x=log10(doses),y=qnval,shape=treatment)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width=0.075) +
  stat_summary(fun.y = mean, geom = "point")+
  facet_wrap(~GF)