## Format Zymo FR D6323 taxonomy table

zymo.GS <- read_delim('./zymo_D6323_data.tsv')

zymo.GS <- as.data.frame(zymo.GS)

count_colname    <- c('Zymo_D6323')

## fix colnames
zymo.GS$lineage <- zymo.GS$Row.names
colnames(zymo.GS)[c(2,3,4,5)] <- c('superkingdom', 'phylum', 'genus', 'species')
##
zymo.GS$count <- zymo.GS[,count_colname]

## drop
zymo.GS <- zymo.GS[,c('superkingdom', 'phylum', 'genus', 'species', 'count', 'lineage')]

## add above taxon to 'unknown'
#add_level <- function(x, level) { x <- paste0(x, '_', level) }
for (level in c('species', 'genus', 'phylum')) {
  #level <- 'phylum'
  i  <- which(colnames(zymo.GS) == level)
  j  <- i - 1
  jc <- colnames(zymo.GS)[j]
  zymo.GS[zymo.GS[,i] == 'unknown', i] <- paste0('unknown_', zymo.GS[zymo.GS[,i] == 'unknown', j])
  #apply(zymo.GS[zymo.GS[,i] == 'unknown', j], 1, function(x) paste0('unknown_', x))
}


fwrite(zymo.GS, './zymo_D6323_data_corrected.tsv')
