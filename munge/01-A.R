# Example preprocessing script.
# get experiment dates ----
experiments <- levels(factor(regmatches(list.files(file.path("data")),
                                        regexpr("[0-9]{6,8}",list.files(file.path("data"))))))


# process data and create datalist ----
datalist.orig <- lapply(experiments, platemeltdown2)
names(datalist.orig) <- experiments

# lets put the list to cache ####
cache('datalist.orig')

# New munge script ----
library(magrittr)
datalist <- "[0-9]{6,8}" %>% 
  Pathfun %>% 
  unique %>% {
    d <- sapply(.,{.%>%paste0("CellGlo_",.,"(_IVIS)?.csv") %>% Pathfun}) %>%
      lapply({.%>% Mungefun %>% melt(id.vars=c("col","row","plate"),variable.name="exp.id")})
    md <- sapply(.,{.%>% paste0("CellGlo_",.,"_[Mm]etadata_([[:alpha:]]+(_GF)?).csv") %>% Pathfun}) %>%
      lapply({.%>% Mungefun}) %>% 
      lapply({.%>%{
        cnam <- colnames(.)%>%
          sub("[A-z]+_[0-9]{6}_[Mm]etadata_([A-z]+(_GF)?).csv","\\1",.)
        date <- colnames(.)%>%
          sub("[A-z]+_([0-9]{6})_[Mm]etadata_([A-z]+(_GF)?).csv","\\1",.)%>%
          as.numeric%>%
          na.omit%>%
          unique
        set_colnames(.,cnam) %>% data.frame(date,.)
      }})
    mapply(inner_join,md,d)
  }