## @knitr fig4 ---

library('ProjectTemplate')
rm(list=ls())
load.project()

# query list object colnames for GF to find out "Sync" experiments
df.sirna <- datalist.orig[names(datalist.orig)%in%c("141008","141022")] %>%
  lapply(., rbindlist) %>% 
  rbindlist

# add instrument types, add row and column names ----
df <- df.sirna %>% 
  mutate(Instrument=ifelse(grepl("IVIS",exp.id),"IVIS","Tecan"),
         rowname=sub("([A-P]{1})([0-9]*)","\\1",well),
         colname=sub("([A-P]{1})([0-9]*)","\\2",well),
         treatment=factor(treatment,levels=c("UT","SCR","siVIM","siCD44")),
         GF=factor(GF,levels=c("bFGF","VEGF", "GDF-2", "FBS_20","FBS_5"), 
                   labels=c("bFGF","VEGF", "GDF-2", "20%FBS","5%FBS")))




library(scales)
df <- df %>% 
  group_by(exp.id,GF) %>% 
  mutate(value = rescale(value))

zfactor <- function(data,treated,untreated){
  t <- data$value[data$treatment%in%treated]
  ut <- data$value[data$treatment%in%untreated]
  1-((3*sd(ut)+3*sd(t))/abs(mean(t)-mean(ut)))
}


# unique(df$GF)

# df %>%
#   group_by(GF,Instrument) %>%
#   do(zfactor(.,treated=c("siCD44"),untreated=c("UT","siSCR","siVIM")))

df %>%
  ggplot(aes(x=value)) + 
  geom_histogram() +
  facet_grid(GF~exp.id,scales="free")



  
library(broom)
mod <- df %>%
  filter(Instrument=="Tecan") %>%
  group_by(GF,doses_GF) %>%
  do(lm(value~treatment,data=.)%>%aov%>%tidy) %>%
  data.frame

# mod <- df %>%
#   filter(Instrument=="Tecan") %>%
#   mutate(induction=paste0(GF,"_",doses_GF)) %>%
#   group_by(induction) %>%
#   do(tukey=lm(value~treatment,data=.)%>%aov%>%TukeyHSD%>%"[["(1))
# 
# situkey <- mod$tukey %>% 
#   lapply(.%>%{data.frame(treatment=rownames(.),.)}) %>% 
#   set_names(mod$induction) %>% 
#   ldply %>%
#   filter(p.adj<=0.05)

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
    ylab(bquote(atop(Relative~cell~number,Mean %+-% SD))) +
    expand_limits(y = 0)
  }


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

# plotlist <- df %>% 
#   dlply(., "GF", function(x) plotlittlelesstuff.mod(x, "IVIS"))
# legend<-g_legend(plotlist[[1]])
# p <- arrangeGrob(arrangeGrob(arrangeGrob(plotlist$bFGF + theme(legend.position="none"),
#                  plotlist$VEGF + theme(legend.position="none"),
#                  plotlist$"20%FBS",
#                  plotlist$"5%FBS"),
#                  legend, ncol=2, widths=c(8,1)),
#                  arrangeGrob(rectGrob(gp=gpar(col="white")),
#                              plotlist$"GDF-2",
#                              rectGrob(gp=gpar(col="white")), 
#                              ncol=3, widths=c(1,4,1)),
#                  nrow=2, heights=c(2,1))
# p
# ggsave(file = paste0("graphs/Sync_siRNA_raw_IVIS_", Sys.Date(),".png"), p, width = 8)


plotlist <- df %>% 
  dlply(., "GF", function(x) plotlittlelesstuff.mod(x, "Tecan"))
legend<-g_legend(plotlist[[1]])
p <- arrangeGrob(arrangeGrob(arrangeGrob(plotlist$bFGF + theme(legend.position="none"),
                                         plotlist$VEGF + theme(legend.position="none"),
                                         plotlist$"20%FBS",
                                         plotlist$"5%FBS") ,
                             legend, ncol=2, widths=c(8,2)),
                 arrangeGrob(rectGrob(gp=gpar(col="white")),
                             plotlist$"GDF-2",
                             rectGrob(gp=gpar(col="white")), 
                             ncol=3, widths=c(1,4,1)),
                 nrow=2, heights=c(2,1))
p
# ggsave(file = paste0("graphs/Sync_siRNA_raw_Tecan_", Sys.Date(),".png"), p, width = 6)
