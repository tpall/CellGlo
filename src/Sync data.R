library('ProjectTemplate')
rm(list=ls())
load.project()
library(magrittr)
# query list object colnames for GF to find out "Sync" experiments
snc <- datalist.orig %>% lapply(., function(x) lapply(x, function(y) sum("GF" %in% names(y)))) %>%
  lapply(., function(z) sum(unlist(z)))
datalist.snc <- datalist.orig[snc > 0]
datalist.snc <- datalist.snc[!names(datalist.snc)%in%c("141008","141022")] # remove siRNA experiments 
names(datalist.snc)

# load data from list and subtract baseline using helper "poolallstuff" ----
df <- datalist.snc %>% poolallstuff

filter(df, is.na(doses))

# log transform values
df$doses <- (df$doses + 0.1)/1e9
df$Instrument <-  "Tecan"
df$Instrument[grep("IVIS", df$exp.id)] <- "IVIS"  
# # add row and column names ----
df$rowname <- sub("([A-P]{1})([0-9]*)","\\1", df$well)
df$colname <- sub("([A-P]{1})([0-9]*)","\\2", df$well)

##### siia peax tulema filtri asi #######
# filter using ivis plate image ----

w96 <- read.csv("data/ivis_plate_filter/140908_IVIS_96w_tuhiplaat.csv", header=FALSE) %>% 
  as.matrix() %>% c() %>% scale(., center=FALSE)

w384 <- read.csv("data/ivis_plate_filter/140908_IVIS_384w_tuhiplaat.csv", header=FALSE) %>% 
  as.matrix() %>% c() %>% scale(., center=FALSE)

df[exp.id == "140610_IVIS",]$value <- df[exp.id == "140610_IVIS",]$value/w96
df[exp.id == "140903_IVIS",]$value <- df[exp.id == "140903_IVIS",]$value/c(w96,w96)
index <- c(outer(LETTERS[1:16], seq(1:24), FUN=paste0)) %in% c(outer(LETTERS[1:16], seq(7,18), FUN=paste0))
df[exp.id == "140705_IVIS",]$value <- df[exp.id == "140705_IVIS",]$value/w384[index,]

df <- df %>% filter(!treatment=="media")

# boxplots from each experiment -- are readouts similar  ----
library(ggplot2)
# intensities in response to different growth factor stimulation in ivis and tecan data
qplot(x = Instrument, y = value, data=df, geom = "boxplot", fill = GF) + 
  facet_wrap(~Instrument, scales = "free")

qplot(x = colname, y = value, data=df, geom = "boxplot", fill = GF) + 
  facet_wrap(~Instrument, scales = "free")

qplot(x = rowname, y = value, data=df, geom = "boxplot", fill = GF) + 
  facet_wrap(~Instrument, scales = "free")


# remove between experiment batch differences ----
# ## Not a good choice actually!!!
# df <- ddply(df, c("Instrument", "GF"), transform, value_ds = detrend(value, exp.id))

# scale/normalise data ----
df <- ddply(df,.(exp.id,GF),mutate,norm.value=scale(value))

# boxplot(cbind(df$norm.value_ds, df$norm.value))

qplot(x = exp.id, y = value, data=df, geom = "boxplot", fill = GF) + 
  facet_wrap(~Instrument, scales = "free")
qplot(x = exp.id, y = norm.value, data=df, geom = "boxplot", fill = GF) + 
  facet_wrap(~Instrument, scales = "free")
# qplot(x = exp.id, y = value_ds, data=df, geom = "boxplot", fill = GF) + 
#   facet_wrap(~Instrument, scales = "free")
# qplot(x = exp.id, y = norm.value_ds, data=df, geom = "boxplot", fill = GF) + 
#   facet_wrap(~Instrument, scales = "free")

grouping <- .(doses,treatment,Instrument,GF,doses_GF)
summary <- ddply(df, grouping, summarise,
                 Mean = mean(value),
                 SD = sd(value),
                 N = length(value),
                 SE = SD/sqrt(N))

summary$treatment2 <- with(summary, paste0(GF, "-", doses_GF, "+\n", treatment))

esteetika <- aes(x = log10(doses), y = Mean, ymin = Mean - SE, ymax = Mean + SE, colour = treatment2)
p <- ggplot(summary, esteetika)
p +  geom_point(size = 3, stat = "identity") + 
     geom_line(size = 1) +
     geom_errorbar(width=0.2) + 
     facet_grid(Instrument~GF, scales = "free") +
     ylab(expression(Mean %+-% SE)) +
     scale_colour_discrete(name = "Treatment")
ggsave(file=paste0("graphs/Snc_cellgrowth_exp_summary_raw_", Sys.Date(),".pdf"))

head(df)

# apply L2 filter to remove plate edge effect row-wise ----
# df <- ddply(df, c("Instrument"), transform, 
#                    filt.value = l2filter.sparse(value, 2))
# 
# df <- ddply(df, c("rowname", "Instrument"), transform, 
#                    normvalue = value/mean(filt.value, na.rm=TRUE))
# 
# summary2 <- ddply(df, grouping, summarize,
#                   Mean = mean(normvalue),
#                   SD = sd(normvalue),
#                   N = length(normvalue),
#                   SE = SD/sqrt(N))
# 
# summary2$treatment2 <- with(summary2, paste0(GF, "-", doses_GF, "+\n", treatment))
# 
# esteetika <- aes(x = log10(doses), y = Mean, ymin = Mean - SE, ymax = Mean + SE, colour = treatment2)
# p <- ggplot(summary2, esteetika)
# p +  geom_point(size = 3, stat = "identity") + 
#   geom_line(size = 1) +
#   geom_errorbar(width=0.2) + 
#   facet_grid(Instrument~GF, scales = "free") +
#   ylab(expression(Mean %+-% SE)) +
#   scale_colour_discrete(name = "Treatment")
# ggsave(paste0("graphs/Sync_cellgrowth_norm", Sys.Date(),".pdf"))


###############################
# summary3 <- ddply(df, grouping, summarize,
#                   Mean = mean(normvalue),
#                   SD = sd(normvalue),
#                   N = length(normvalue),
#                   SE = SD/sqrt(N))
# 
# summary3$treatment2 <- with(summary3, paste0(GF, "-", doses_GF, "+\n", treatment))


# plotstuff <- function(input, error){
#   library(ggthemes)
#   esteetika <- eval(parse(text = paste("aes(x = log10(doses), y = Mean, ymin = Mean -", 
#                                        error,
#                                        ", ymax = Mean +",
#                                        error, 
#                                        ", colour = treatment2)")))
#   ggplot(input, esteetika) +
#     geom_point(size = 3, stat = "identity") + 
#     geom_line(size = 1) +
#     geom_errorbar(width=0.2) + 
#     ylab(parse(text=paste("Mean %+-%", error))) +
#     facet_wrap(~Instrument, scales = "free") +
#     scale_colour_colorblind(name = "Treatment")
# }

# plotlist <- dlply(summary3, "GF", function(x) plotstuff(x, error = "SD"))
# names(plotlist)
# g <- arrangeGrob(plotlist$bFGF,plotlist$VEGF,plotlist$"GDF-2",plotlist$HGF)
# g
# ggsave(file = paste0("graphs/Sync_cellgrowth_norm_breakup", Sys.Date(),".pdf"), g, width = 8)

# plotlittlelesstuff <- function(input, error, instr){
#   library(dplyr)
#   library(ggthemes)
#   esteetika <- eval(parse(text = paste("aes(x = log10(doses), y = Mean, ymin = Mean -", 
#                                        error,
#                                        ", ymax = Mean +",
#                                        error, 
#                                        ", colour = treatment2)")))
#   ggplot(filter(input, Instrument == instr), esteetika) +
#     geom_point(size = 3, stat = "identity") + 
#     geom_line(size = 1) +
#     geom_errorbar(width=0.2) + 
#     ylab(parse(text=paste("Mean %+-%", error))) +
#     scale_colour_colorblind(name = "Treatment")
# }

# plotlist <- dlply(summary3, "GF", function(x) plotlittlelesstuff(x, error = "SD"))
# g <- arrangeGrob(plotlist$bFGF,plotlist$VEGF,plotlist$"GDF-2",plotlist$HGF)
# ggsave(file = paste0("graphs/Sync_cellgrowth_norm_breakup_Tecan_", Sys.Date(),".pdf"), g, width = 8)
# g

# NB! uses de-trended data and introduces noise!
summary4 <- ddply(df, grouping, summarize,
                  Mean = mean(value_ds),
                  SD = sd(value_ds),
                  N = length(value_ds),
                  SE = SD/sqrt(N))

summary4$treatment2 <- with(summary4, paste0(GF, "-", doses_GF, "+\n", treatment))

plotlist <- dlply(summary4, "GF", function(x) plotlittlelesstuff(x, "SD", "IVIS"))
g <- arrangeGrob(plotlist$bFGF,plotlist$VEGF,plotlist$"GDF-2",plotlist$HGF)
ggsave(file = paste0("graphs/Sync_cellgrowth_raw_breakup_IVIS_", Sys.Date(),".pdf"), g, width = 8)
g

summary4 <- ddply(df, grouping, summarize,
                  Mean = mean(norm.value_ds),
                  SD = sd(norm.value_ds),
                  N = length(norm.value_ds),
                  SE = SD/sqrt(N))


summary4$treatment2 <- with(summary4, paste0(GF, "-", doses_GF, "+\n", treatment))
plotlist <- summary4 %>% 
  filter(doses<1.60001e-05) %>% 
  dlply(., "GF", function(x) plotlittlelesstuff(x, "SD", "IVIS"))
g <- arrangeGrob(plotlist$bFGF,plotlist$VEGF,plotlist$"GDF-2",plotlist$HGF)
ggsave(file = paste0("graphs/Sync_cellgrowth_norm.value_ds_IVIS_", Sys.Date(),".pdf"), g, width = 8)
g

plotlist <- summary4 %>% 
  filter(doses<1.60001e-05) %>% 
  dlply(., "GF", function(x) plotlittlelesstuff(x, "SD", "Tecan"))
g <- arrangeGrob(plotlist$bFGF,plotlist$VEGF,plotlist$"GDF-2",plotlist$HGF)
ggsave(file = paste0("graphs/Sync_cellgrowth_norm.value_ds_Tecan_", Sys.Date(),".pdf"), g, width = 8)
g

plotlist <- summary4 %>% 
  filter(doses<1.60001e-05) %>% 
  dlply(., "GF", function(x) plotstuff(x, "SD"))
g <- arrangeGrob(plotlist$bFGF,plotlist$VEGF,plotlist$"GDF-2",plotlist$HGF)
ggsave(file = paste0("graphs/Sync_cellgrowth_norm.value_ds_", Sys.Date(),".pdf"), g, width = 8)
g




# qplot(y = rowname, x = as.numeric(colname), 
#       data = ddply(df, "Instrument", mutate, value = scale(value)), 
#       size = value, geom = "point", alpha = value) + 
#       facet_wrap(~Instrument, scales = "free")
# ggsave("graphs/plateimage.pdf")

# NB! onl scaled data, no 
grouping <- c("doses", "treatment", "Instrument" , "GF", "doses_GF")
summary <- ddply(df, grouping, summarize,
                  Mean = mean(norm.value),
                  SD = sd(norm.value),
                  N = length(norm.value),
                  SE = SD/sqrt(N))

summary$treatment2 <- with(summary, paste0(GF, "-", doses_GF, "+\n", treatment))

plotlist <- summary %>% 
  filter(doses<1.60001e-05) %>% 
  dlply(., "GF", function(x) plotlittlelesstuff(x, "SD", "IVIS"))
g <- arrangeGrob(plotlist$bFGF,plotlist$VEGF,plotlist$"GDF-2",plotlist$HGF)
ggsave(file = paste0("graphs/Sync_cellgrowth_norm.value_IVIS_", Sys.Date(),".pdf"), g, width = 8)
g

plotlist <- summary %>% 
  filter(doses<1.60001e-05) %>% 
  dlply(., "GF", function(x) plotlittlelesstuff(x, "SD", "Tecan"))
g <- arrangeGrob(plotlist$bFGF,plotlist$VEGF,plotlist$"GDF-2",plotlist$HGF)
ggsave(file = paste0("graphs/Sync_cellgrowth_norm.value_Tecan_", Sys.Date(),".pdf"), g, width = 8)
g

plotlist <- summary %>% 
  filter(doses<1.60001e-05) %>% 
  dlply(., "GF", function(x) plotstuff(x, "SD"))
g <- arrangeGrob(plotlist$bFGF,plotlist$VEGF,plotlist$"GDF-2",plotlist$HGF)
ggsave(file = paste0("graphs/Sync_cellgrowth_norm.value_", Sys.Date(),".pdf"), g, width = 8)
g


library(ggthemes)
df$treatment2 <- with(df, paste0(GF, "-", doses_GF, "+\n", treatment))
df %>%
  filter(doses<1.60001e-05) %>%
  filter(Instrument == "Tecan") %>%
  filter(GF == "HGF") %>%
  ggplot(., aes(x=log10(doses), y=norm.value, color = treatment2)) +
  geom_point() +
  geom_smooth() + 
  scale_colour_colorblind(name = "Treatment")

# GDF-2 data ---

library(ggthemes)
df$treatment2 <- with(df, paste0(GF, "-", doses_GF, "+\n", treatment))
gdf <- df %>% 
  filter(doses<1.60001e-05) %>%
  filter(Instrument == "Tecan") %>%
  mutate(value=scale(value))

p <- gdf %>%
  filter(GF == "GDF-2"&treatment%in%c("rhIgG-Fc","3MUT-Fc")) %>%
  ggplot(., aes(x=log10(doses), y = value, color = treatment2)) +
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x), 
               fun.ymax = function(x) mean(x) + sd(x), 
               geom = "pointrange", size=1) +
  stat_summary(fun.y = mean, geom = "line") +
  scale_colour_colorblind(name = "Treatment") +
  xlab("Protein concentration, Log(M)") + ylab(expression(Mean %+-% SD))
q <- gdf %>%
  filter(GF == "GDF-2"&treatment%in%c("UT")) %>%
  ggplot(., aes(x=factor(doses_GF), y = value)) +
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x), 
               fun.ymax = function(x) mean(x) + sd(x), 
               geom = "pointrange", size=1) +
  #stat_summary(fun.y = mean, geom = "bar") +
  xlab("GDF-2, ng/ml") + 
  ylab(NULL) + ggtitle("GDF-2 is growth\nsuppressive to EC")
r<-arrangeGrob(q,p,nrow=1, widths=c(1, 4))
r
r<-p + annotation_custom(ggplotGrob(q), xmin = -6.75, xmax = -5.6, ymin=-2.25, ymax=-1.5)
ggsave("graphs/GDF2_stimulated_HUVEC_3MUT-Fc.png", r, width=6, dpi=300)
