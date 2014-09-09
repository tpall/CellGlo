library('ProjectTemplate')
load.project()

# query list object colnames for GF to find out "Sync" experiments
snc <- lapply(lapply(datalist.orig, function(x) lapply(x, function(y) sum("GF" %in% names(y)))), function(z) sum(unlist(z)))
datalist.snc <- datalist.orig[snc > 0]
names(datalist.snc)

# load data from list and subtract baseline using helper "poolallstuff" ----
df <- ldply(datalist.snc, poolallstuff)

# log transform values
df$doses <- (df$doses + 0.1)/1e9
df$Instrument <-  "Tecan"
df$Instrument[grep("IVIS", df$exp.id)] <- "IVIS"  
# # add row and column names ----
df$rowname <- sub("([A-P]{1})([0-9]*)","\\1", df$well)
df$colname <- sub("([A-P]{1})([0-9]*)","\\2", df$well)



# boxplots from each experiment -- are readouts similar  ----
library(ggplot2)
qplot(x = Instrument, y = value, data=df, geom = "boxplot", fill = GF) + 
  facet_wrap(~Instrument, scales = "free")

qplot(x = colname, y = value, data=df, geom = "boxplot", fill = GF) + 
  facet_wrap(~Instrument, scales = "free")

qplot(x = rowname, y = value, data=df, geom = "boxplot", fill = GF) + 
  facet_wrap(~Instrument, scales = "free")

# scale/normalise data ----
df <- ddply(df, c("exp.id", "GF"), mutate, norm.value = scale(value, center = FALSE))
boxplot(df$norm.value)

qplot(x = exp.id, y = norm.value, data=df, geom = "boxplot", fill = GF) + 
  facet_wrap(~Instrument, scales = "free")

grouping <- c("doses", "treatment", "Instrument" , "celldensity", "GF", "doses_GF")
summary <- ddply(df, grouping, summarize,
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
df <- ddply(df, c("Instrument"), transform, 
                   filt.value = l2filter.sparse(value, 2))

df <- ddply(df, c("rowname", "Instrument"), transform, 
                   normvalue = value/mean(filt.value, na.rm=TRUE))

summary2 <- ddply(df, grouping, summarize,
                  Mean = mean(normvalue),
                  SD = sd(normvalue),
                  N = length(normvalue),
                  SE = SD/sqrt(N))

summary2$treatment2 <- with(summary2, paste0(GF, "-", doses_GF, "+\n", treatment))

esteetika <- aes(x = log10(doses), y = Mean, ymin = Mean - SE, ymax = Mean + SE, colour = treatment2)
p <- ggplot(summary2, esteetika)
p +  geom_point(size = 3, stat = "identity") + 
  geom_line(size = 1) +
  geom_errorbar(width=0.2) + 
  facet_grid(Instrument~GF, scales = "free") +
  ylab(expression(Mean %+-% SE)) +
  scale_colour_discrete(name = "Treatment")
ggsave(paste0("graphs/Sync_cellgrowth_norm", Sys.Date(),".pdf"))


###############################
summary3 <- ddply(df, grouping, summarize,
                  Mean = mean(norm.value),
                  SD = sd(norm.value),
                  N = length(norm.value),
                  SE = SD/sqrt(N))

summary3$treatment2 <- with(summary3, paste0(GF, "-", doses_GF, "+\n", treatment))

plotstuff <- function(input){
  library(ggthemes)
  esteetika <- aes(x = log10(doses), y = Mean, ymin = Mean - SE, ymax = Mean + SE, colour = treatment2)
  ggplot(input, esteetika) +
    geom_point(size = 3, stat = "identity") + 
    geom_line(size = 1) +
    geom_errorbar(width=0.2) + 
    facet_wrap(~Instrument, scales = "free") +
    ylab(expression(Mean %+-% SE)) +
    scale_colour_colorblind(name = "Treatment")
}

plotlist <- dlply(summary3, "GF", plotstuff)
names(plotlist)
g <- arrangeGrob(plotlist$bFGF,plotlist$VEGF,plotlist$"GDF-2",plotlist$HGF)
ggsave(file = paste0("graphs/Sync_cellgrowth_norm_breakup", Sys.Date(),".pdf"), g, width = 8)

plotlittlelesstuff <- function(input){
  library(dplyr)
  library(ggthemes)
  esteetika <- aes(x = log10(doses), y = Mean, ymin = Mean - SE, ymax = Mean + SE, colour = treatment2)
  ggplot(filter(input, Instrument == "Tecan"), esteetika) +
    geom_point(size = 3, stat = "identity") + 
    geom_line(size = 1) +
    geom_errorbar(width=0.2) + 
    ylab(expression(Mean %+-% SE)) +
    scale_colour_colorblind(name = "Treatment")
}

plotlist <- dlply(summary3, "GF", plotlittlelesstuff)
g <- arrangeGrob(plotlist$bFGF,plotlist$VEGF,plotlist$"GDF-2",plotlist$HGF)
ggsave(file = paste0("graphs/Sync_cellgrowth_norm_breakup_Tecan_", Sys.Date(),".pdf"), g, width = 8)


qplot(y = rowname, x = as.numeric(colname), 
      data = ddply(df, "Instrument", mutate, value = scale(value)), 
      size = value, geom = "point", alpha = value) + 
      facet_wrap(~Instrument, scales = "free")
ggsave("graphs/plateimage.pdf")
