# improved function ----
platemeltdown2 <- function(experiment){
  require(plyr)
  # read file list
  datafile <- regmatches(list.files(file.path("data")),
                         regexpr(paste0("CellGlo_", experiment,
                                        "(_IVIS)?.csv"), list.files(file.path("data"))))
  # x <- datafile[1]
  # create data list
  make.datalist <- function(x){
    value <- read.csv(file.path("data", x), header=FALSE)
    data.frame(value = c(as.matrix(value)))
  }
    
  realdata <- alply(datafile, 1, make.datalist)
  names(realdata) <- sub(paste0("CellGlo_(", experiment, "(_IVIS)?).csv"), "\\1", datafile)
  
  # metadata magic
  metadatafiles <- regmatches(list.files(file.path("data")), regexpr(
    paste0("CellGlo_", experiment,
           "_[Mm]etadata_([[:alpha:]]+(_GF)?).csv"), list.files(file.path("data"))))
  
  metadatanames <- sub(paste0("CellGlo_", experiment, "_[Mm]etadata_([[:alpha:]]+(_GF)?).csv"), "\\1", metadatafiles)
  
  extractmetadata <- function(x){
    metadata <- read.csv(x, header=FALSE) # metadata files don't have header, just matrix 
    wells <- c(outer(LETTERS[1:nrow(metadata)], seq(ncol(metadata)), FUN=paste , sep=""))
    data.frame(well = wells, value = c(as.matrix(metadata))[1:length(wells)])
  }
  
  # merge data and metadata
  m <- lapply(file.path("data", metadatafiles), extractmetadata)
  m <- join_all(m, "well") # plyr function for joining files in list
  colnames(m) <- c("well", metadatanames)
  
  # add metadata to values
  realdata <- lapply(realdata, function(x) data.frame(m, value = x))
  mapply(function(x, y) {data.frame(x, exp.id = as.factor(y))}, realdata, names(realdata), SIMPLIFY = FALSE)
  }

### New munge functions ####
library(magrittr)
library(plyr);library(dplyr)
Pathfun <- . %>% 
  regexpr(.,list.files(file.path("data"))) %>%
  regmatches(list.files(file.path("data")),.)

Mungefun <- . %>%
  sapply({.%>%file.path("data",.)%>%
            read.csv(header=FALSE)%>%
            as.matrix%>%t%>%c})%>%{
              di <- dim(.) %>% "["(1)
              if(di==384) id <- expand.grid(col=c(1:24),row=LETTERS[1:16],plate=seq(1))
              if(di%in%c(96,192)) id <- expand.grid(col=c(1:12),row=LETTERS[1:8],plate=seq(di/96))
              data.frame(id,.)}
