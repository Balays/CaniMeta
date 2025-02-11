
#library(data.table)

pattern  <- '_EMU_EMUdb'
res.dirs <- grep(pattern, list.dirs('.', T, F), value = T)

res.dirs <- c("./FR_D6323_EMU_EMUdb", "./MCM_D6300_EMU_EMUdb", "./MCM_D6331_EMU_EMUdb")

res.dir  <- res.dirs[1] #'MCM_D6331_Emu_EMUdb'

plot_all <- data.table(NULL)

for (res.dir in res.dirs[]) {
  try({
  file     <- paste0(res.dir, '/Detection statistics.tsv')
  message('Importing: ', file)
  plot_sum <- fread(file = file)
  plot_sum[,source := file]
  #plot_sum   <- rbind(correl.exclunclass, correl.keepunclass)
  print(plot_sum[,.N,sample])

  plot_all <- rbind(plot_sum, plot_all, fill=T)
  })
}


plot_all[,level := factor(level, levels = c('Genus', 'Species'))]
plot_all[,GS_version := factor(GS_version, levels=c('D6300', 'D6331', 'D6323'))]

plot_all[,.N,statistic]
plot_all[,.N,GS_version]


### Plotting

res.dir <- "MCM_Detect_stats"
dir.create(res.dir)


## plot Chi2-statistics

plot_chi <- data.table(pivot_wider(
  plot_all[level == 'Species' & (statistic == 'Chi_stat' | statistic == 'Chisq.test.sign'), ],
  names_from = 'statistic', values_from = 'value'
))


ggChi <- ggplot(plot_chi[,] ,
                 aes(x = sample_name, y = Chi_stat, fill = primer)) +
  geom_col(position = position_dodge2()) +
  geom_text(data=plot_chi[Chisq.test.sign == 1,],
            aes(x = sample_name, y = Chi_stat + 10, fill = primer),
            label = '*', position = position_dodge2()) +
  guides(fill = guide_legend(title = "Source", direction = "horizontal", nrow = 1)) +
  scale_fill_manual(values=as.character(pal.man[])) +  # c(10:16)
  scale_y_continuous(name='Chi-statistic')+
  #facet_wrap(~ sample, scales = "free") +
  #coord_cartesian(ylim = c(0,1)) +
  theme_ipsum() +
  theme(legend.position="bottom",
        axis.text.x = element_text(angle=-60),
        axis.title.x = element_blank(),
        panel.spacing.x = unit(0.1, 'mm')) +
  facet_nested(cols=vars(GS_version),
               rows=vars(Detection_threshold_percent),
               # ncol=5,
               scales = 'free') +
  ggtitle('Chi2 statistics by detection threshold')


ggsave(file.path(res.dir, 'Chi-stat-combined.jpg'), ggChi, height = 16, width = 16)




## plot 'F1', 'Precision' and 'Recall' Statistics

plot_sum <- plot_all[statistic %in% c('F1', 'Precision', 'Recall'),
                     .(mean_value=mean(value), sd_value=sd(value)),
                     by=.(Detection_threshold_percent, statistic, GS_version, primer, level)]

ggall <- ggplot(plot_sum,
                aes(x = GS_version, y = mean_value, fill = primer)) +
  geom_col(position = position_dodge2()) +
  #geom_text(#data=plot_chi[Chisq.test.sign == 1,],
  #          aes(x = sample_name, y = statistic, fill = primer),
  #          label = Chisq.test.sign, position = position_dodge2()) +
  geom_errorbar(aes(x=GS_version, group=primer, ymin=mean_value - sd_value, ymax=mean_value + sd_value),
                position = position_dodge2()) +
  guides(fill = guide_legend(title = "Source", direction = "horizontal", nrow = 1)) +
  scale_fill_manual(values=as.character(pal.man[])) +  # c(10:16)
  #scale_y_continuous(name='Chi-statistic')+
  #facet_wrap(~ sample, scales = "free") +
  #coord_cartesian(ylim = c(0,1)) +
  theme_ipsum() +
  theme(legend.position="bottom",
        axis.text.x = element_text(angle=-60),
        axis.title.x = element_blank(),
        panel.spacing.x = unit(0.1, 'mm')) +
  facet_nested(cols=vars(Detection_threshold_percent, level),
               rows=vars(statistic),
               # ncol=5,
               scales = 'free') +
  ggtitle('Detection statistics')

ggsave(file.path(res.dir, 'Detect-stats-combined.jpg'), ggall, height = 16, width = 20)



