library(drc)
fit <- tec %>%
  filter(!treatment=="media") %>%
  filter(doses<1.50001e-05) %>%
  group_by(date,plate,GF) %>% 
  mutate(value=(value-min(value))/(range(value)%>%diff)) %>%
  filter(treatment%in%c("rhIgG-Fc","3MUT-Fc")) %>%
  mutate(value=value/median(value[treatment=="rhIgG-Fc"])) %>%
  group_by(date,GF,doses,treatment,treat2) %>%
  summarise(value=mean(value)) %>%
  filter(GF=="bFGF") %>%
  drm(value~doses,treat2,data=.,
      fct=LL.4(names = c("b", "lower", "upper", "ed50")),
      control=drmc(method = "Nelder-Mead"))

tec %>%
  filter(!treatment=="media") %>%
  filter(doses<1.50001e-05) %>%
  group_by(date,plate,GF) %>% 
  mutate(value=(value-min(value))/(range(value)%>%diff)) %>%
  filter(treatment%in%c("rhIgG-Fc","3MUT-Fc")) %>%
  mutate(value=value/median(value[treatment=="rhIgG-Fc"])) %>%
  group_by(date,GF,doses,treatment,treat2) %>%
  summarise(value=mean(value)) %>%
  filter(GF=="bFGF") %>%
  ggplot(aes(x=log10(doses),y=value,shape=treat2)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width=0.075) +
  stat_summary(fun.y = mean, geom = "point") +
  geom_line(aes(y=predict(fit))) +
  ylim(0,1.2)



Mytidy <- . %>% summary %>% coef %>% {
  term <- rownames(.)
  data.frame(term,., row.names = NULL) %>% 
    set_colnames(c("term","estimate","std.error","t.value","p.value"))
}

fit%>%Mytidy
predict(fit)
coef(fit) %>% names

models <- tec %>%
  filter(!treatment=="media") %>%
  filter(doses<1.50001e-05) %>%
  group_by(date,plate,GF) %>% 
  mutate(value=(value-min(value))/(range(value)%>%diff)) %>%
  filter(treatment%in%c("rhIgG-Fc","3MUT-Fc")) %>%
  mutate(value=value/median(value[treatment=="rhIgG-Fc"])) %>%
  group_by(date,GF,doses,treatment,treat2) %>%
  summarise(value=mean(value)) %>%
  filter(GF=="HGF") %>%
  bootstrap(1000) %>%
  do(mod=drm(value~doses,curveid=treat2,data=.,
               fct=LL.4(names = c("b", "lower", "upper", "ed50")),
      control=drmc(method = "Nelder-Mead",errorm=FALSE,maxIt=1000,rmNA=FALSE)))

bootcoef <- models %>% do(data.frame(coefs=coef(.$mod)))
cfs <- c("b:25 bFGF 3MUT-Fc","b:25 bFGF rhIgG-Fc","lower:25 bFGF 3MUT-Fc",
         "lower:25 bFGF rhIgG-Fc","upper:25 bFGF 3MUT-Fc","upper:25 bFGF rhIgG-Fc",
         "ed50:25 bFGF 3MUT-Fc","ed50:25 bFGF rhIgG-Fc")
bootcoef %>% { 
  rep <- nrow(.) %>% divide_by(8) %>% rep(cfs,.)
  data.frame(rep,.)
} %>% 
  group_by(rep) %>%
  
  
  summarise(Mean=mean(coefs,na.rm = TRUE),
            lower=quantile(coefs,0.025,na.rm = TRUE),
            upper=quantile(coefs,0.975,na.rm = TRUE))
