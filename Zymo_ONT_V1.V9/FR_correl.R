
## check the correlation between the experimntal results for FR and the experimental results of Hsapiens and Clupus


library(readr)
library(data.table)

zymo.GS <- read_delim('../standards/zymo_6300_composition.txt')


FR_GS <- fread('FR_D6323_EMU_EMUdb/FR_D6323_EMU_EMUdb_all_species_ratios.tsv')
FR_GS <- FR_GS[,-c(1:6,8)]

FR_GS <- melt(FR_GS, id.vars = 1, variable.name = 'sample', value.name = 'count')
