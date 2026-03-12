library(data.table)
library(stringi)

metadata  <- fread('PRJEB82097_meta.tsv')
metadata[, sample := biosample]
setnames(metadata, c('Sample_name', 'Sex', 'Date'), c('name', 'sex', 'date'))
metadata <- data.frame(metadata, row.names = metadata$sample)

emu.dir   <- 'EMU/rel_abundance/min1500.max1800'

emu.files <- list.files(emu.dir, 'rel-abundance.tsv', full.names = T)

file <- emu.files[1]

read_emu <- function(file, sample_pattern = '.lenfilt.*') {

  message('importing ', file, '...')

  dt <- fread(file, na.strings = '')

  sample <- basename(file)
  sample <- gsub(sample_pattern, '', sample)

  dt[, lineage := apply(.SD, 1, function(x) {
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x)) paste(x, collapse = "::") else NA_character_
  }), .SDcols = ranks]
  dt[is.na(lineage), lineage := tax_id]

  dt[, sample := sample]

  return(dt)

}

# ranks

make_phylo <- function(emu_dt,
                       metadata,
                       ranks = c("superkingdom", "phylum", "class", "order", "family", "genus", "species"),
                       unassigned = c('mapped_filtered', 'mapped_unclassified', 'unmapped') ) {

  # Set these to 'unassigned'
  if( !is.null(unassigned) ) {

    emu_unassigned <- emu_dt[lineage %in% unassigned, ]
    emu_unassigned[,lineage := 'unassigned']
    emu_unassigned <- emu_unassigned[,.(`estimated counts` = sum(`estimated counts`)), by = .(sample, lineage)]

    emu_dt  <- emu_dt[!lineage %in% unassigned, ]
    emu_dt  <- rbind(emu_dt, emu_unassigned, fill = T)
  }
  cols   <- c('lineage', ranks)
  taxtab <- unique(emu_dt[,..cols])
  #taxtab[,"tax_id" := NULL]
  taxtab <- data.frame(taxtab[,], row.names = taxtab$lineage)


  otutab <- dcast(emu_dt, lineage ~ sample, value.var = 'estimated counts', fill = 0)
  otutab <- data.frame(otutab[,-1], row.names = otutab$lineage)
  #colnames(otutab) <- sample

  #otutab <- otutab[rownames(taxtab), ]

  #sampdat <- data.frame( unique(emu_all_dt[,.(sample)]) ) ;  row.names(sampdat) <- sampdat$sample

  metadata <- as.data.frame(metadata)
  sampdat  <- metadata[metadata$sample %in% unique(emu_dt$sample), ]

  ps <- phyloseq(tax_table(as.matrix(taxtab)), otu_table(as.matrix(otutab),taxa_are_rows = T), sample_data(sampdat))

  return(ps)

}


## All data
emu_all_dt <- rbindlist(purrr::map(emu.files, read_emu))
ps.all     <- make_phylo(emu_all_dt, metadata)



## Save PS Objects
saveRDS(ps.all , 'EMU/ps/min1500.max1800/serteperti.exp.all.PS.rds')






