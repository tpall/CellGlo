library('ProjectTemplate')
load.project()

mitmes <- 1
# load data
df <- datalist.orig[[mitmes]][[1]]
exp.date <- names(datalist.orig)[mitmes]
# remove baseline
baseline <- df[df$treatment=="media",]
df <- df[!df$treatment=="media",]
df$value <- df$value - mean(baseline$value)

# log transform values
df$doses <- (df$doses+0.1)/1e9

# filter out 2e+5
# df <- df[!df$value < 2e+5,]

p <- ggplot(df, aes(x = factor(signif(doses, 2)), y = value, fill = treatment)) +
  geom_boxplot() + facet_wrap(~time) + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

summary <- ddply(df, c("doses", "treatment", "time"), summarize,
      Mean = mean(value),
      SD = sd(value),
      N = length(value))

q <- ggplot(summary, aes(x = log10(doses), y = Mean, colour = treatment)) +
  geom_point(size = 3) + facet_wrap(~time) + geom_line(size = 1) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width=0.2) 

library(gridExtra)
g <- arrangeGrob(p, q, ncol=1)
ggsave(file=paste0("graphs/cellgrowth_", exp.date,".pdf"), g)

# lets check some column or row wise trends ----
# row labels are in "rowname" and column labels are in "colname"
rows <- qplot(x = rowname, y = value, geom = 'boxplot', data = df) + facet_wrap(~time)
cols <- qplot(x = colname, y = value, geom = 'boxplot', data = df) + facet_wrap(~time, scales = "free_x")

g <- arrangeGrob(rows, cols, ncol=1)
ggsave(file=paste0("graphs/platetrend_", exp.date,".pdf"), g)

# L2 filtering approach ----
library(Matrix)
# http://www.exegetic.biz/blog/2014/03/filtering-data-with-l2-regularisation/
#   ?utm_source=rss&utm_medium=rss&utm_campaign=filtering-data-with-l2-regularisation
# copy paste L2 filter for sparse matrixes
l2filter.sparse <- function(x, lambda = 0.0) {
  n <- length(x)
  I = Diagonal(n)
  D = bandSparse(n = n - 2, m = n, k = c(0, 1, 2),
                 diagonals = list(rep(1, n), rep(-2, n), rep(1, n)))
  (solve(I + 2 * lambda * t(D) %*% D) %*% x)[,1]
}

# apply L2 filter to remove plate edge effect row-wise
df <- ddply(df, "time", transform, filt.value = scale(l2filter.sparse(value, 2), center = FALSE))
ggplot(df, aes(x = colname, y = rowname, size = value, colour = filt.value)) + geom_point()
r <- qplot(x = rowname, y = filt.value, geom = 'boxplot', data = df) + facet_wrap(~time) + ggtitle("Raw data")
s <- qplot(x = colname, y = filt.value, geom = 'boxplot', data = df) + facet_wrap(~time, scales = 'free_x') + ggtitle("L2 filter values")
g <- arrangeGrob(r, s, ncol=1)
ggsave(file=paste0("graphs/smoothing_", exp.date,".pdf"), g)

df <- ddply(df, c("time", "rowname"), transform, normvalue = value/mean(filt.value))

# df$normvalue <- scale(df$normvalue, center=FALSE)
qplot(x = rowname, y = normvalue, geom = 'boxplot', data = df) + facet_wrap(~time)
qplot(x = colname, y = normvalue, geom = 'boxplot', data = df) + facet_wrap(~time, scales = 'free_x')

# calculate and remove bottom (bottom == 24 hour timepoints) ----
bottom <- mean(df$normvalue[df$time==24])
df$normvalue <- df$normvalue - bottom
# iqrnorm <- function(x){
#   m <- median(x, na.rm = T)
#   iqr <- IQR(x, na.rm = T)
#   (x - m)/iqr
# }
df <- ddply(df, c("time"), transform, normvalue = scale(normvalue, center = FALSE))

summary2 <- ddply(df, c("doses","treatment", "time"), summarize,
                 Mean = mean(normvalue),
                 SD = sd(normvalue),
                 N = length(normvalue))

ggplot(summary2[summary2$time==72,], aes(x = log10(doses), y = Mean, colour = treatment)) +
  geom_point(size = 3) + facet_wrap(~time) + geom_line(size = 1) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width=0.2) 
ggsave(paste0("graphs/cellgrowth_norm", exp.date,".pdf"))

ggplot(df[df$time==72,], aes(x = factor(signif(doses, 2)), y = normvalue, colour = treatment)) +
  geom_point(size = 3) + stat_smooth(aes(group = treatment), method="lm", size = 1)
ggsave("graphs/filtered_data_lm-smooth.pdf")

# fit function to data ----
fits <- dlply(df, c("time","treatment"), function(x) lm(normvalue~doses, data = x))
lapply(fits, summary)
lapply(fits, anova)

library(drc)
fits <- dlply(df[df$time==72,], "treatment", function(x) drm(normvalue~doses, data = x, fct = LL.5()))

# andmed tuleb min-max normaliseerida!!!!
plot(fits[[1]])
plot(fits[[3]])
plot(fits[[2]])
dev.off()
