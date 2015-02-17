helper.function <- function()
{
  return(1)
}

# function to load data from list and subtract baseline ----
poolallstuff <- function(datalist){
  library(data.table)
  library(magrittr)
  library(dplyr)
  # load data and in case we used also IVIS to quantitate we add extra ldply step
  datalist %>% lapply(., rbindlist) %>% rbindlist %>%
    group_by(exp.id) %>% 
    mutate(value = value - mean(value[treatment == "media"], na.rm = TRUE))
    #filter(., !treatment=="media")
}


# good old detrend function
detrend <- function(var, exp){resid(lm(var~factor(exp))) + mean(var, na.rm = TRUE)}

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

## plotstuff----
plotstuff <- function(input, error){
  library(ggthemes)
  esteetika <- eval(parse(text = paste("aes(x = log10(doses), y = Mean, ymin = Mean -", 
                                       error,
                                       ", ymax = Mean +",
                                       error, 
                                       ", colour = treatment2)")))
  ggplot(input, esteetika) +
    geom_point(size = 3, stat = "identity") + 
    geom_line(size = 1) +
    geom_errorbar(width=0.2) + 
    ylab(parse(text=paste("Mean %+-%", error))) +
    facet_wrap(~Instrument, scales = "free") +
    scale_colour_colorblind(name = "Treatment")
}

#plotlittlelesstuff----
plotlittlelesstuff <- function(input, error, instr){
  library(dplyr)
  library(ggthemes)
  esteetika <- eval(parse(text = paste("aes(x = log10(doses), y = Mean, ymin = Mean -", 
                                       error,
                                       ", ymax = Mean +",
                                       error, 
                                       ", colour = treatment2)")))
  ggplot(filter(input, Instrument == instr), esteetika) +
    geom_point(size = 3, stat = "identity") + 
    geom_line(size = 1) +
    geom_errorbar(width=0.2) + 
    ylab(parse(text=paste("Mean %+-%", error))) +
    scale_colour_colorblind(name = "Treatment")
}