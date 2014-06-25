helper.function <- function()
{
  return(1)
}

# function to load data from list and subtract baseline ----
poolallstuff <- function(datalist){
  # load data and in case we used also IVIS to quantitate we add extra ldply step
  df <- ldply(datalist)
  # remove baseline
  df <- ddply(df, c("exp.id"), transform, value = value - mean(value[treatment == "media"]))
  df[!df$treatment=="media", ]
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