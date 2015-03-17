library('ProjectTemplate')
rm(list=ls())
load.project()
library(magrittr)
# query list object colnames for GF to find out "Sync" experiments
datalist.sirna <- datalist.orig[names(datalist.orig)%in%c("141008","141022")]
names(datalist.sirna)

# load data from list and subtract baseline using helper "poolallstuff" ----
df <- datalist.sirna %>% lapply(., rbindlist) %>% rbindlist %>% 
  group_by(exp.id) %>% mutate(value = scale(value, center=FALSE))

df$Instrument <-  "Tecan"
df$Instrument[grep("IVIS", df$exp.id)] <- "IVIS"  
# # add row and column names ----
df$rowname <- sub("([A-P]{1})([0-9]*)","\\1", df$well)
df$colname <- sub("([A-P]{1})([0-9]*)","\\2", df$well)
df$treatment <- factor(df$treatment , levels=c("UT","SCR","siVIM","siCD44"))
df$GF <- factor(df$GF , levels=c("bFGF","VEGF", "GDF-2", "FBS_20","FBS_5"), 
                labels=c("bFGF","VEGF", "GDF-2", "20%FBS","5%FBS"))
unique(df$GF)

#plotlittlelesstuff----
plotlittlelesstuff.mod <- function(input, instr){
  library(dplyr)
  library(ggthemes)
  library(Hmisc)
  gfname <- unique(input$GF)
  if(length(unique(input$doses_GF))==1){
    setup <- aes(x=treatment,y=value, color=treatment)
  } else {setup <- aes(x=factor(doses_GF), y=value, color=treatment, group=treatment)}
  
  p <- ggplot(filter(input, Instrument == instr), setup)
  
  if(length(unique(input$doses_GF))==1){
    p <- p + stat_summary(fun.ymin = function(x) mean(x) - sd(x), 
                          fun.ymax = function(x) mean(x) + sd(x), 
                          geom = "errorbar",width=0.25) +
      stat_summary(fun.y = mean,
                   geom = "point") +
      scale_colour_colorblind(guide = FALSE)+
      xlab("Treatment")
  } else {
    p <- p + stat_summary(fun.y = mean,
                          fun.ymin = function(x) mean(x) - sd(x), 
                          fun.ymax = function(x) mean(x) + sd(x), 
                          geom = "pointrange") +
      stat_summary(fun.y = mean,
                   geom = "line") +
      scale_colour_colorblind(name = "Treatment") +
      xlab("Growth factor\nconcentration")
  }
  
  p  + ggtitle(gfname) + 
    ylab(expression(Mean %+-% SD))
  
}
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

plotlist <- df %>% 
  dlply(., "GF", function(x) plotlittlelesstuff.mod(x, "IVIS"))
legend<-g_legend(plotlist[[1]])
p <- arrangeGrob(arrangeGrob(arrangeGrob(plotlist$bFGF + theme(legend.position="none"),
                 plotlist$VEGF + theme(legend.position="none"),
                 plotlist$"20%FBS",
                 plotlist$"5%FBS"),
                 legend, ncol=2, widths=c(8,1)),
                 arrangeGrob(rectGrob(gp=gpar(col="white")),
                             plotlist$"GDF-2",
                             rectGrob(gp=gpar(col="white")), 
                             ncol=3, widths=c(1,4,1)),
                 nrow=2, heights=c(2,1))
p
ggsave(file = paste0("graphs/Sync_siRNA_raw_IVIS_", Sys.Date(),".png"), p, width = 8)


plotlist <- df %>% 
  dlply(., "GF", function(x) plotlittlelesstuff.mod(x, "Tecan"))
p <- arrangeGrob(arrangeGrob(arrangeGrob(plotlist$bFGF + scale_y_continuous(limits=c(0.5,1.75)) + theme(legend.position="none"),
                                         plotlist$VEGF + scale_y_continuous(limits=c(0.5,1.75)) + theme(legend.position="none"),
                                         plotlist$"20%FBS" + scale_y_continuous(limits=c(0.5,1.75)),
                                         plotlist$"5%FBS"  + scale_y_continuous(limits=c(0.5,1.75))) ,
                             legend, ncol=2, widths=c(8,2)),
                 arrangeGrob(rectGrob(gp=gpar(col="white")),
                             plotlist$"GDF-2" + scale_y_continuous(limits=c(0.5,1.75)),
                             rectGrob(gp=gpar(col="white")), 
                             ncol=3, widths=c(1,4,1)),
                 nrow=2, heights=c(2,1))
p
ggsave(file = paste0("graphs/Sync_siRNA_raw_Tecan_", Sys.Date(),".png"), p, width = 6)
