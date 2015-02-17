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

###########################
# function to read data into dataframe and add metadata, produces datalist ----
# experiment <- experiments[6]
# platemeltdown <- function(experiment){
#   require(plyr)
#   # read file list
#   datafile <- regmatches(list.files(file.path("data")),
#                          regexpr(paste0("CellGlo_", experiment,
#                                         ".csv"), list.files(file.path("data"))))
#   
#   z <- read.csv(file.path("data", datafile), header=FALSE)
#   colnames(z) <- seq(ncol(z))
#   rownames(z) <- paste(LETTERS[1:nrow(z)])
#   # wells <- c(outer(LETTERS[1:nrow(z)], seq(ncol(z)), FUN=paste , sep=""))
#   z$rowname <- rownames(z)
#   z <- melt(z, id.vars = "rowname", variable_name = "colname")
#   z$well <- paste0(z$rowname, z$colname)
#   
#   # metadata magic
#   metadatafiles <- regmatches(list.files(file.path("data")), regexpr(
#     paste0("CellGlo_", experiment,
#            "_[Mm]etadata_([[:alpha:]]+)_*([[:alpha:]]*).csv"), list.files(file.path("data"))))
#   
#   extractmetadatanames <- function(x) {
#     m <- regexec("_[Mm]etadata_([[:alpha:]]+)_*([[:alpha:]]*)",x)
#     parts <- do.call(rbind,
#                      lapply(regmatches(x, m), `[`, c(2L,3L)))
#     parts <- apply(parts, 1, paste, collapse = ".")
#     parts <- gsub("\\.$","",parts)
#     parts <- c(na.omit(parts))
#     parts}
#   
#   metadatanames <- extractmetadatanames(metadatafiles)
#   
#   extractmetadata <- function(x, wells){
#     md <- read.csv(x, header=FALSE) # metadata files don't have header, just matrix 
#     metadata <- data.frame(well = wells, value = c(as.matrix(md))[1:length(wells)])
#     metadata}
#   
#   # merge data and metadata
#   m <- lapply(file.path("data", metadatafiles), extractmetadata, wells = z$well)
#   m <- join_all(m, "well") # plyr function
#   colnames(m) <- c("well", metadatanames)
#   z <- join(m, z)
#   z$exp.id <- as.factor(experiment)
# 
# z}
# experiment <- experiments[6]