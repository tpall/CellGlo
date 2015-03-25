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

# query list object col numbers to find "Sync" experiments;snc=12,async=10,sirna=9 ----
snc <- datalist %>% 
  lapply({.%>%names%>%length}) %>% 
  unlist %>% 
  equals(12) %>% 
  datalist[.] %>%
  bind_rows

snc %<>% mutate(doses = doses%>%as.numeric,
                doses_GF = doses_GF%>%as.numeric,
                GF = GF%>%as.factor,
                date = date%>%as.factor)

# log transform treatment doses ----
snc %<>% mutate(doses = (doses+0.1)/1e9)

# lets use only Tecan measurements ----
tec <- snc %>% filter(!grepl("IVIS",exp.id))
ivi <- snc %>% filter(grepl("IVIS",exp.id))

# add variables content for normalisation and treat2
tec %<>%
  mutate(content = ifelse(treatment=="UT","pos","sample"),
         content = ifelse(doses_GF==0,"neg",content),
         content = ifelse(treatment=="media","blank",content),
         content = ifelse(treatment=="FUM","xFUM",content),
         treat2 = paste(doses_GF,GF,treatment))
