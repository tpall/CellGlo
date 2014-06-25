df <- datalist.orig[[2]]

# remove baseline
baseline <- df[df$treatment=="media",]
df <- df[!df$treatment=="media",]
df$value <- df$value - median(baseline$value)

# filter out 2e+5
df <- df[!df$value < 2e+5,]

p <- ggplot(df, aes(x = factor(doses), y = value, fill = treatment)) +
  geom_boxplot() + facet_wrap(~time) + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

summary <- ddply(df, c("doses","treatment", "time"), summarize,
                 Mean = mean(value),
                 SD = sd(value),
                 N = length(value))

q <- ggplot(summary, aes(x = log(doses+0.1), y = Mean, colour = treatment)) +
  geom_point(size = 3) + geom_line(size = 1) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width=0.2) + facet_wrap(~time)

library(gridExtra)
g <- arrangeGrob(p, q, ncol=1)
g
ggsave(file="graphs/cellgrowth_ivis.pdf", g)


