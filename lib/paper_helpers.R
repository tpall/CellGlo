# plot summary data
Plotfun <- .%>% {
  ggplot(.,aes(x=log10(doses),y=value,color=treatment)) +
    stat_summary(fun.data = mean_se, geom = "pointrange") +
    stat_summary(fun.y = mean, geom = "line",size=1)+
    scale_color_colorblind()}

# edge effect removal functions
Predfun <- function(x) x$col[2:3]%>%
  data.frame(col=c(2,3),eff=.)%>%
  lm(eff~col,data=.)%>%
  predict(.,data.frame(col=1))

Medfun <- function(x,newfirst) {
  x$col[1]<-newfirst
  x}