lapply(datalist.snc, function(x) sapply(x, colnames))

lapply(datalist.snc, function(x) lapply(x, function(y) apply(y, 2, unique)))

unique(summary2$treatment2)
