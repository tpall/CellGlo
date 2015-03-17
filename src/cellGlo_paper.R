library('ProjectTemplate')
rm(list=ls())
load.project()

# query list object colnames for GF to find out "Sync" experiments
snc <- datalist %>% 
  lapply(., function(x) lapply(x, function(y) sum("GF" %in% names(y)))) %>%
  lapply(., function(z) sum(unlist(z))) %>% { snc <- unlist(.) %>% is_greater_than(0)
                                              datalist.orig[snc]} 
snc <- snc[!names(snc)%in%c("141008","141022")] # remove siRNA experiments 

# Extract Tecan data
tec <- snc %>% lapply("[[",1)
tec %>% lapply(summary)

tec %>% lapply({.%>%mutate(value=scale(value,center = FALSE)%>%c)%>%
                  filter(treatment=="media")})


filter(df, is.na(doses))

# log transform values
df$doses <- (df$doses + 0.1)/1e9
df$Instrument <-  "Tecan"
df$Instrument[grep("IVIS", df$exp.id)] <- "IVIS"  
# # add row and column names ----
df$rowname <- sub("([A-P]{1})([0-9]*)","\\1", df$well)
df$colname <- sub("([A-P]{1})([0-9]*)","\\2", df$well)