# Example preprocessing script.
# get experiment dates ----
experiments <- levels(factor(regmatches(list.files(file.path("data")),
                                        regexpr("[0-9]{6,8}",list.files(file.path("data"))))))


# process data and create datalist ----
datalist.orig <- lapply(experiments, platemeltdown2)
names(datalist.orig) <- experiments

# lets put the list to cache ####
cache('datalist.orig')
