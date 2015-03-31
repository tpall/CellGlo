library('ProjectTemplate')
load.project()

# query list object colnames for GF to find out "async" experiments
asnc <- lapply(lapply(datalist.orig, function(x) lapply(x, function(y) sum(!"GF" %in% names(y)))), function(z) sum(unlist(z)))
datalist.asnc <- datalist.orig[asnc > 0]

# load data from list and subtract baseline using helper "poolallstuff" ----
df <- ldply(datalist.asnc, rbind_all)
df %<>% filter(!grepl("IVIS",exp.id))

# log transform values
df$doses <- (df$doses+0.1)/1e9
# df$Instrument <-  "Tecan"
# df$Instrument[grep("IVIS", df$exp.id)] <- "IVIS"  
# # add row and column names ----
df$rowname <- sub("([A-P]{1})([0-9]*)","\\1", df$well)
df$colname <- sub("([A-P]{1})([0-9]*)","\\2", df$well)

# in 140509 experiment IVIS data were aquired only at 72 hour timepoint
# qplot(x = factor(time), y = value, data = df[df$exp.id=="140509_IVIS",], 
# geom = "boxplot")
# lets remove 24 and 48 hour data from this experiment
# df <- df[!(df$exp.id == "140509_IVIS" & df$time %in% c("24", "48")), ]
df %<>% filter(time==72)

# boxplots from each experiment -- are readouts similar  ----
qplot(x = factor(exp.id), y = value, data=df, geom = "boxplot") 

# ok, let's check if different experiments values differ
# fits <- dlply(threedays, c("Instrument"), 
#               function(x) lm(value~factor(exp.id), data = x))
# lapply(fits, anova)
df %>% 
  lm(value~factor(exp.id), data = .) %>%
  anova

df %>%
  group_by(exp.id) %>%
  mutate(value=scale(value)%>%c) %T>%
  qplot(x = factor(exp.id), y = value, data=., geom = "boxplot") %>% 
  lm(value~factor(exp.id), data = .) %>%
  anova

df %>%
  group_by(exp.id) %>%
  mutate(value=scale(value)%>%c) %>%
  mutate(value=(value-min(value))/(range(value)%>%diff)) %>%
    ggplot(aes(x=log10(doses),y=value,shape=treatment)) +
    stat_summary(fun.data = mean_se, geom = "pointrange") +
    stat_summary(fun.y = mean, geom = c("line","point"))

# apply L2 filter to remove plate edge effect row-wise ----
df %>% filter(treatment=="FUM") %>% select(doses) %>% unique
df %>% filter(treatment=="rhIgG-Fc") %>% select(doses) %>% unique

df %>%
  filter(!exp.id=="140512") %>%
  group_by(exp.id) %>% 
  filter(!treatment=="media") %>%
  group_by(exp.id) %>% 
  mutate(value = (value-min(value))/(range(value)%>%diff)) %>% {
    fumctrl <- filter(.,treatment=="rhIgG-Fc"&doses%in%c(1.28000e-08,4.01000e-08)) %>%
      mutate(doses=2.01e-08)
    rbind(.,fumctrl)
  } %>%
  group_by(exp.id,doses) %>%
  mutate(value = value/median(value[treatment=="rhIgG-Fc"], na.rm = TRUE)) %>%
  filter(!(treatment=="rhIgG-Fc"&doses==2.01e-08)) %>%
  filter(!(doses==1.60001e-05|treatment=="UT")) %>%
  mutate(treatment = factor(treatment,c("rhIgG-Fc","3MUT-Fc","FUM"))) %>%
  ggplot(aes(x = log10(doses), y = value, shape = treatment)) +
  stat_summary(fun.data = mean_se, geom = "pointrange") +
  stat_summary(fun.y = mean, geom = "point", size=3) +
  stat_summary(fun.y = mean, geom = "line") +
  scale_y_continuous(expand = c(0, 0)) +
  expand_limits(y = 0) + 
  theme_classic() +
  theme(legend.justification=c(1,0), legend.position=c(1,0), legend.title=element_blank()) +
  ylab("Relative cell number") +
  xlab(bquote(list(Conc.,log[10](M))))
