



#krona.files <- list.files('sourmash_NCBI', '.k31.kreport.txt', full.names = T)

sampdat.all = data.frame(#platform=c('MySeq', 'NovaSeq', 'ONT'),
                         Vregion='WGS',
                         workflow='sourmash',
                         method='k31',
                         db='NCBI')

ps  <- import.sourmash.files(indir = 'Toti_sourmash_NCBI', sampdat.all = sampdat.all)

taxtab  <- data.frame(ps@tax_table)
otutab  <- data.frame(ps@otu_table)
sampdat <- data.frame(ps@sam_data)


sample_names(ps) <- gsub('_unmapped', '', sample_names(ps))
ps@sam_data$sample <- sample_names(ps)

sampdat <- data.frame(ps@sam_data)

## DNA Isolation Kit
ps@sam_data$DNA.Isolation.Kit <- fifelse(grepl('_I', sample_names(ps)), 'I',
                                         fifelse(grepl('_MN', sample_names(ps)), 'MN',
                                                 fifelse(grepl('_Q', sample_names(ps)), 'Q',
                                                         fifelse(grepl('_Z', sample_names(ps)), 'ZHMW',
                                                                 fifelse(grepl('NovaSeq', sample_names(ps)), 'ZMB',
                                                                         fifelse(grepl('ONT', sample_names(ps)), 'ZMB','NA'))))))

## Sample name
ps@sam_data$sample_name <- stri_extract_all_regex(ps@sam_data$sample, '[Q|Z|I][0-9]')
ps@sam_data$sample_name[is.na(ps@sam_data$sample_name)] <- stri_extract_all_regex(ps@sam_data$sample[is.na(ps@sam_data$sample_name)], '.[MN][0-9]')
ps@sam_data$sample_name <- gsub(' ', '', ps@sam_data$sample_name)
ps@sam_data$sample_name[grepl('NovaSeq|ONT', ps@sam_data$sample)] <- 'ZMB1'

## Sample number
ps@sam_data$sample_nr <- as.integer(gsub('[Q|Z|I|M|N]', '', ps@sam_data$sample_name))
ps@sam_data$sample_nr[grepl('NovaSeq|ONT', ps@sam_data$sample)] <- 1

## Sample ID
sample_names(ps)   <- gsub('_S[0-9]*.*', '', sample_names(ps))
sample_names(ps)   <- gsub('NovaSeq_WGS', 'NovaSeq', sample_names(ps))
ps@sam_data$sample <- sample_names(ps)

## Platform
ps@sam_data$platform <- fifelse(grepl('NovaSeq', ps@sam_data$sample), 'NovaSeq',
                                fifelse(grepl('ONT', ps@sam_data$sample), 'ONT', 'MySeq'))

## OK
taxtab  <- data.frame(ps@tax_table)
otutab  <- data.frame(ps@otu_table)
sampdat <- data.frame(ps@sam_data)


## merge different sequencings of the same sample

combine_ps <- function(ps, merge_vars = c('Vregion', 'workflow', 'db', 'platform', 'sample_name') # , 'DNA.Isolation.Kit', 'sample_nr',
) {

  taxtab  <- data.frame(ps@tax_table)
  otutab  <- data.frame(ps@otu_table, taxon_id=taxa_names(ps)); setDT(otutab)
  sampdat <- data.frame(ps@sam_data); setDT(sampdat); sampdat[,sample := sample_names(ps)]

  ps.data <- melt.data.table(otutab, id.vars='taxon_id', value.name = 'count', variable.name = 'sample')
  ps.data <- merge(sampdat, ps.data, by='sample')

  ## summarise
  #merge_vars  <- c(merge_vars, 'taxon_id')
  ps.data.sum <- ps.data[,.(count = sum(count)), by=c(merge_vars, 'taxon_id')]
  ps.data.sum <- unite(ps.data.sum, 'sample', merge_vars, remove = F, na.rm = T)

  ## spread back
  #ps.data.sp <- spread(ps.data.sum, 'taxon_id', value = 'count')
  ps.data.sp <- spread(ps.data.sum[,c('sample', "taxon_id", 'count')], 'sample', value = 'count')

  ##
  otutab  <- data.frame(ps.data.sp[,-1], row.names = ps.data.sp$taxon_id)

  #taxtab  <- data.frame(ps@tax_table)

  sampdat <- unique(sampdat[,..merge_vars])
  sampdat <- unite(sampdat, 'sample', merge_vars, remove = F, na.rm = T)
  rownames(sampdat) <- sampdat$sample

  otutab <- otutab[,sampdat$sample]
  otutab <- otutab[rownames(taxtab), ]

  ps.merged <- phyloseq(otu_table(as.matrix(otutab), taxa_are_rows = T),
                        tax_table(as.matrix(taxtab)),
                        sample_data(sampdat))

  return(ps.merged)
}

ps.sourmash <- combine_ps(ps, merge_vars = c('Vregion', 'workflow', 'db', 'platform', 'sample_name'))

ps <- ps.sourmash

## OK
taxtab  <- data.frame(ps@tax_table)
otutab  <- data.frame(ps@otu_table)
sampdat <- data.frame(ps@sam_data)


## DNA Isolation Kit
ps@sam_data$DNA.Isolation.Kit <- fifelse(grepl('NovaSeq', sample_names(ps)), 'ZMB',
                                         fifelse(grepl('ONT', sample_names(ps)), 'ZMB',
                                                 fifelse(grepl('_I', sample_names(ps)), 'I',
                                                         fifelse(grepl('_MN', sample_names(ps)), 'MN',
                                                                 fifelse(grepl('_Q', sample_names(ps)), 'Q',
                                                                         fifelse(grepl('_Z', sample_names(ps)), 'ZHMW', 'NA'))))))

## Sample number
ps@sam_data$sample_nr <- as.integer(gsub('[Q|Z|I|M|N]', '', ps@sam_data$sample_name))
ps@sam_data$sample_nr[grepl('NovaSeq|ONT', ps@sam_data$sample)] <- 1

## OK
taxtab  <- data.frame(ps@tax_table)
otutab  <- data.frame(ps@otu_table)
sampdat <- data.frame(ps@sam_data)

ps.sourmash <- ps

saveRDS(ps.sourmash, 'Toti_Illumina_WGS_sourmash.k31_NCBI_PS.rds')







