# lets scale by RMS each plate at each timepoint
mask <- ddply(df, c("exp.id", "time"), transform, value = (value-min(value))/(max(value)-min(value)))
mask <- ddply(mask, "well", summarize, well.mask = mean(value, na.rm = TRUE))
mask$row <- gsub("([A-P]{1})([0-9]{1,2})", "\\1", mask$well)
mask$col <- gsub("([A-P]{1})([0-9]{1,2})", "\\2", mask$well)
mask$row <- factor(mask$row, levels = unique(as.character(mask$row)))
mask$col <- as.numeric(as.character(mask$col))
qplot(x = col, y = row, data = mask, size = well.mask, geom = "point")

for(i in seq_along(mask$well)){
  replace <- which(threedays$well == mask$well[i])
  threedays$norm2[replace] <- threedays$value_db[replace] / mask$well.mask[i]
}

summary <- ddply(threedays, c("doses", "treatment", "Instrument" , "celldensity"), summarize,
                 Mean = mean(norm2),
                 SD = sd(norm2),
                 N = length(norm2),
                 SE = SD/sqrt(N))

q <- ggplot(summary, aes(x = log10(doses), y = Mean, colour = treatment)) +
  geom_point(size = 3) + facet_grid(Instrument~celldensity, scales = "free_y") + geom_line(size = 1) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width=0.2) 
ggsave(file=paste0("graphs/cellgrowth_exp_summary_raw_", Sys.Date(),".pdf"), q)

#
mask.row <- ddply(mask, "row", summarize, well.mask = mean(well.mask))

for(i in seq_along(mask.row$row)){
  replace <- which(threedays$row == mask.row$row[i])
  threedays$norm2[replace] <- threedays$value_db[replace] / mask.row$well.mask[i]
}

summary <- ddply(threedays, c("doses", "treatment", "Instrument" , "celldensity"), summarize,
                 Mean = mean(norm2),
                 SD = sd(norm2),
                 N = length(norm2),
                 SE = SD/sqrt(N))

q <- ggplot(summary, aes(x = log10(doses), y = Mean, colour = treatment)) +
  geom_point(size = 3) + facet_grid(Instrument~celldensity, scales = "free_y") + geom_line(size = 1) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width=0.2) 
ggsave(file=paste0("graphs/cellgrowth_exp_summary_raw_", Sys.Date(),".pdf"), q)