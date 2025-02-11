

pattern  <- '_EMU_EMUdb'
res.dirs <- grep(pattern, list.dirs('.', T, F), value = T)

res.dirs <- c("./FR_D6323_EMU_EMUdb", "./MCM_D6300_EMU_EMUdb", "./MCM_D6331_EMU_EMUdb")

res.dir  <- res.dirs[1] #'MCM_D6331_Emu_EMUdb'

plot_all <- data.table(NULL)


for (res.dir in res.dirs[]) {
  try({
    file     <- paste0(res.dir, '/Observed-Theoretical Correlations.tsv')
    message('Importing: ', file)
    plot_sum <- fread(file = file)
    plot_sum[,source := file]
    #plot_sum   <- rbind(correl.exclunclass, correl.keepunclass)
    print(plot_sum[,.N,sample])

    plot_all <- rbind(plot_sum, plot_all, fill=T)
  })
}

plot_all[,rank := factor(rank, levels = ranks)]

plot_all[,GS_version := factor(GS_version, levels=c('D6300', 'D6331', 'D6323'))]



plot_all  <- unique(plot_all[,.(GS_version, primer, rank, sample_name, R2, source)])
plot_sum  <- plot_all[,
                     .(mean_value=mean(R2), sd_value=sd(R2)),
                     by=.(rank, GS_version, primer)]



### Plot R2 values

res.dir <- "MCM_Detect_stats"
dir.create(res.dir)



### Barplot
ggr2 <- ggplot(plot_sum[rank != 'phylum'],
               aes(x = GS_version, y = mean_value, fill = primer)) +
  geom_col(position = position_dodge2()) +
  geom_errorbar(aes(x=GS_version, group=primer, ymin=mean_value - sd_value, ymax=mean_value + sd_value),
                position = position_dodge2()) +
  guides(fill = guide_legend(title = "Primer", direction = "horizontal", nrow = 1)) +
  scale_fill_manual(values=as.character(pal.man[])) +  # c(10:16)
  #facet_wrap(~ sample, scales = "free") +
  coord_cartesian(ylim = c(0,1)) +
  theme_ipsum() +
  theme(legend.position="bottom",
        axis.text.x = element_text(angle=-60),
        panel.spacing.x = unit(0.1, 'mm')) +
  facet_nested(cols=vars(GS_version),
               rows=vars(rank),
               # ncol=5,
               scales = 'free') +
  ggtitle('R2 values')


ggr2
ggsave(file.path(res.dir, 'R2-combined_barplot.jpg'), ggr2, height = 14, width = 12)


### Lineplot

plot_all[,group := gsub('[0-9]*', '', sample_name)]


ggr2 <- ggplot(plot_all, aes(x=rank, y=R2, fill=sample_name )) +
  geom_line(aes(group = sample_name)) +
  geom_point(shape=21, size=4) +
  guides(fill = guide_legend(direction = "vertical", ncol = 1)) +
  scale_fill_manual(values=as.character(pal.man)) +
  labs(title = "Lineplot of Pearson's Correlation Across Taxon Levels",
       y = "R-squared value",
       x = "Taxon Level",
       fill = "Sample Name") +
  theme_ipsum() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.75, vjust = 0.25)) +
  facet_nested(rows=vars(GS_version), cols=vars(group), scales = 'fixed')


ggr2
ggsave(file.path(res.dir, 'R2-combined_lineplot.jpg'), ggr2, height = 14, width = 12)
