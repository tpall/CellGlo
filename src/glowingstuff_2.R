library('ProjectTemplate')
load.project()

# query list object colnames for GF to find out "async" experiments
asnc <- lapply(lapply(datalist.orig, function(x) lapply(x, function(y) sum(!"GF" %in% names(y)))), function(z) sum(unlist(z)))
datalist.asnc <- datalist.orig[asnc > 0]

# load data from list and subtract baseline using helper "poolallstuff" ----
df <- ldply(datalist.asnc, poolallstuff)

# log transform values
df$doses <- (df$doses+0.1)/1e9
df$Instrument <-  "Tecan"
df$Instrument[grep("IVIS", df$exp.id)] <- "IVIS"  
# # add row and column names ----
df$rowname <- sub("([A-P]{1})([0-9]*)","\\1", df$well)
df$colname <- sub("([A-P]{1})([0-9]*)","\\2", df$well)

# in 140509 experiment IVIS data were aquired only at 72 hour timepoint
# qplot(x = factor(time), y = value, data = df[df$exp.id=="140509_IVIS",], 
# geom = "boxplot")
# lets remove 24 and 48 hour data from this experiment
df <- df[!(df$exp.id == "140509_IVIS" & df$time %in% c("24", "48")), ]

# boxplots from each experiment -- are readouts similar  ----
qplot(x = factor(exp.id), y = value, data=df, geom = "boxplot") + 
  facet_grid(Instrument~time, scales = "free_y")

# lets work on with 72 hour datapoint ----
threedays <- df[df$time==72,]

# ok, let's check if different experiments values differ
fits <- dlply(threedays, c("Instrument"), 
              function(x) lm(value~factor(exp.id), data = x))
lapply(fits, anova)

# plot data at 72 hour timepoint
qplot(x = factor(exp.id), y = value, data=threedays, geom = "boxplot") + 
  facet_grid(Instrument~time, scales = "free_y")

# lets detrend data then -----
threedays <- ddply(threedays, "Instrument", 
                   transform, value_db = detrend(value, exp.id))

# was detrending effective
fits <- dlply(threedays, c("Instrument"), 
              function(x) lm(value_db ~ factor(exp.id), data = x))
lapply(fits, anova) # yep seems so ...
qplot(x = factor(exp.id), y = value_db, data=threedays, geom = "boxplot") + 
  facet_grid(Instrument~time, scales = "free_y")
# 

summary <- ddply(threedays, c("doses", "treatment", "Instrument" , "celldensity"), 
                 summarize,
                 Mean = mean(value_db),
                 SD = sd(value_db),
                 N = length(value_db),
                 SE = SD/sqrt(N))

q <- ggplot(summary, aes(x = log10(doses), y = Mean), colour = treatment) +
  geom_point(size = 3) + facet_grid(Instrument~celldensity, scales = "free_y") + 
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width=0.2) + 
  ylab(expression(Mean %+-% SE))
ggsave(file=paste0("graphs/Asnc_cellgrowth_exp_summary_raw_", Sys.Date(),".pdf"), q)
# 
p <- ggplot(threedays, aes(x = factor(signif(doses, 2)), y = value_db, fill = treatment)) +
  geom_boxplot() + facet_wrap(~Instrument, scales = "free_y") + 
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

library(gridExtra)
g <- arrangeGrob(p, q, ncol=1)
ggsave(file=paste0("graphs/Asnc_cellgrowth_exp_summary_", Sys.Date(),".pdf"), g)

# tegelikult peaks arvutama reakeskmised ja nende pealt see reafunktsioon ----


# apply L2 filter to remove plate edge effect row-wise ----
threedays <- ddply(threedays, c("exp.id"), transform, 
            filt.value = l2filter.sparse(value_db, 2))

# lets look at data after filtering ----
ggplot(threedays, aes(x = factor(signif(doses, 2)), y = filt.value, fill = treatment)) +
  geom_boxplot() + facet_wrap(~Instrument, scales = "free_y") + 
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

threedays <- ddply(threedays, c("rowname", "exp.id"), transform, 
            normvalue = value_db/mean(filt.value, na.rm=TRUE))

# lets look at data after adjusting ----
ggplot(threedays, aes(x = factor(signif(doses, 2)), y = normvalue, fill = treatment)) +
  geom_boxplot() + facet_wrap(~Instrument, scales = "free_y") + 
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

# df$normvalue <- scale(df$normvalue, center=FALSE)
#qplot(x = rowname, y = normvalue, geom = 'boxplot', data = threedays) + facet_grid(Instrument~time, scales = "free_y")
# qplot(x = colname, y = normvalue, geom = 'boxplot', data = threedays) + facet_grid(Instrument~time, scales = 'free_y')


summary2 <- ddply(threedays, c("doses","treatment", "Instrument"), 
                  summarize,
                  Mean = mean(normvalue),
                  SD = sd(normvalue),
                  N = length(normvalue),
                  SE = SD/sqrt(N))

ggplot(summary2, aes(x = log10(doses), y = Mean, colour = treatment)) +
  geom_point(size = 3) + 
  facet_wrap(~Instrument, scales = "free_y") + 
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width=0.2) + 
  ylab(expression(Mean %+-% SE))
ggsave(paste0("graphs/Asnc_cellgrowth_norm", Sys.Date(),".pdf"))

# ggplot(threedays, aes(x = factor(signif(doses, 2)), y = normvalue, colour = treatment)) +
#   geom_point(size = 3) + stat_smooth(aes(group = treatment), method="lm", size = 1) +
#   facet_wrap(~Instrument)
# ggsave(paste0("graphs/filtered_data_lm-smooth", exp.date,".pdf"))
# 
# # fit function to data ----
# fits <- dlply(df, c("time","treatment", "Instrument"), function(x) lm(normvalue~doses, data = x))
# lapply(fits, summary)
# lapply(fits, anova)