
## Import Libraries

library(tidyverse)
library(ggpubr)
library(factoextra)
library(FactoMineR)
library(ggrepel)
library(phyloseq, quietly = T)
library(DESeq2)
library(data.table)
library(ggsci)
library(hrbrthemes)
library(ggh4x)
library(RColorBrewer)
library(stringr); library(stringi)
library("ape")
library("ggtree")
library(Biostrings)
library(dplyr)
library(UpSetR)
library(ComplexHeatmap)
library(genefilter)
library(viridis)

## Import own functions -->> SET PATH !!!
misc.dir    <- 'C:/GitHub/Rlyeh/R'
minitax.dir <- 'C:/GitHub/minitax/R'

if (.Platform$OS.type!="windows") {

  misc.dir  <- paste0('/mnt/', gsub(':', '/', tolower(gsub('/.*', '', misc.dir))),
                      stri_replace_first_regex(misc.dir, '.*:\\/', '')
  )

  minitax.dir  <- paste0('/mnt/', gsub(':', '/', tolower(gsub('/.*', '', minitax.dir))),
                         stri_replace_first_regex(minitax.dir, '.*:\\/', '')
  )

}

for(f in list.files(misc.dir,    '*.R', full.names = T)) { try({message(f); source(f)}) }
for(f in list.files(minitax.dir, '*.R', full.names = T)) { try({message(f); source(f)}) }

bam.flags           <- fread(paste0(misc.dir, '/bam.flags.tsv'))
gff_compare_classes <- fread(paste0(misc.dir, '/gff_compare.txt'))
nproc <- 32


#### IMPORT KRAKEN2 OUTPUTS #####

#krona.files <- list.files('sourmash_NCBI', '.k31.kreport.txt', full.names = T)

sampdat.all = data.frame(#platform=c('MySeq', 'NovaSeq', 'ONT'),
                         host ='Toti',
                         Vregion='WGS',
                         workflow='kraken2',
                         #method='k31',
                         db='kraken2_viral')

indir <- 'Toti_host_removed_kraken2_viral_db'
ps  <- import.krona.files(indir = indir, krona.files = list.files(indir, pattern = '*_krona.txt', full.names = T, recursive = T), sampdat.all = sampdat.all)

taxtab  <- data.frame(ps@tax_table)
otutab  <- data.frame(ps@otu_table)
sampdat <- data.frame(ps@sam_data)


sample_names(ps) <- gsub('_unmapped', '', sample_names(ps))
ps@sam_data$sample <- sample_names(ps)

sampdat <- data.frame(ps@sam_data)

## DNA Isolation Kit
ps@sam_data$DNA.Isolation.Kit <- fifelse(grepl('NovaSeq', sample_names(ps)), 'ZMB',
                                         fifelse(grepl('ONT', sample_names(ps)), 'ZMB',
                                                 fifelse(grepl('_I', sample_names(ps)), 'I',
                                                         fifelse(grepl('_MN', sample_names(ps)), 'MN',
                                                                 fifelse(grepl('_Q', sample_names(ps)), 'Q',
                                                                         fifelse(grepl('_Z', sample_names(ps)), 'ZHMW', 'NA'))))))


## Sample name
ps@sam_data$sample_name <- stri_extract_all_regex(ps@sam_data$sample, '[Q|Z|I][0-9]')
ps@sam_data$sample_name[is.na(ps@sam_data$sample_name)] <- stri_extract_all_regex(ps@sam_data$sample[is.na(ps@sam_data$sample_name)], '.[MN][0-9]')
ps@sam_data$sample_name <- gsub(' ', '', ps@sam_data$sample_name)
ps@sam_data$sample_name[grepl('NovaSeq|ONT', ps@sam_data$sample)] <- 'Z'

## Sample number
ps@sam_data$sample_nr <- as.integer(gsub('[Q|Z|I|M|N]', '', ps@sam_data$sample_name))
ps@sam_data$sample_nr[grepl('NovaSeq|ONT', ps@sam_data$sample)] <- 1

ps@sam_data$sample_name <- paste0(ps@sam_data$DNA.Isolation.Kit, '', ps@sam_data$sample_nr)

## Sample ID
sample_names(ps)   <- gsub('_S[0-9]*.*', '', sample_names(ps))
sample_names(ps)   <- gsub('NovaSeq_WGS', 'NovaSeq', sample_names(ps))
sample_names(ps)   <- gsub('ONT_WGS', 'ONT', sample_names(ps))
ps@sam_data$sample <- sample_names(ps)

## Platform
ps@sam_data$platform <- fifelse(grepl('NovaSeq', ps@sam_data$sample), 'NovaSeq', fifelse(grepl('ONT', ps@sam_data$sample), 'ONT', 'MySeq'))

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

ps.kraken2 <- combine_ps(ps, merge_vars = c('Vregion', 'workflow', 'db', 'platform', 'sample_name', 'host'))

ps <- ps.kraken2

## DNA Isolation Kit
ps@sam_data$DNA.Isolation.Kit <- fifelse(grepl('NovaSeq', sample_names(ps)), 'ZMB',
                                         fifelse(grepl('ONT', sample_names(ps)), 'ZMB',
                                                 fifelse(grepl('_I', sample_names(ps)), 'I',
                                                         fifelse(grepl('_MN', sample_names(ps)), 'MN',
                                                                 fifelse(grepl('_Q', sample_names(ps)), 'Q',
                                                                         fifelse(grepl('_Z', sample_names(ps)), 'ZHMW', 'NA'))))))

## Sample number
ps@sam_data$sample_nr <- as.integer(gsub('[[:upper:]]*', '', ps@sam_data$sample_name))
#ps@sam_data$sample_nr[grepl('NovaSeq', ps@sam_data$sample)] <- 1


## Sample ID
ps@sam_data$sample <- paste0(ps@sam_data$host, '_', ps@sam_data$platform, '_', ps@sam_data$sample_name)
sample_names(ps)   <- ps@sam_data$sample


## OK
taxtab  <- data.frame(ps@tax_table)
otutab  <- data.frame(ps@otu_table)
sampdat <- data.frame(ps@sam_data)

ps.kraken2 <- ps


#### SAVE KRAKEN2 OUTPUTS #####

saveRDS(ps.kraken2, 'Toti_Illumina_WGS_kraken2_Viral_PS.rds')







