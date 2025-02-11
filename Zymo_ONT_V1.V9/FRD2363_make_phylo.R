# Load necessary libraries
library(phyloseq)
library(data.table)
library(stringr)

# Read in the Kraken2 output
kraken_data <- fread("../standards/Zymo.D6323.DNA.report/in2165.DNA/ZymoBIOMICS.fecal.ref.illumina.pe/AbundanceTables/AbundanceTable.csv")

# Fix column name issue if necessary
setnames(kraken_data, "#OTU ID", "OTU_ID")

# Parse taxonomy and abundance data
tax_data <- kraken_data[, tstrsplit(OTU_ID, split = ";", fixed = TRUE)]
setnames(tax_data, c("Kingdom", "Phylum", "Genus", "Species", "X", "Y", "Z"))

# Check taxon levels to ensure completeness
all_taxa_levels <- c("Kingdom", "Phylum", "Genus", "Species", "X", "Y", "Z")
missing_taxa <- setdiff(all_taxa_levels, names(tax_data))
if (length(missing_taxa) > 0) {
  stop(paste("The following taxon levels are missing:", paste(missing_taxa, collapse = ", ")))
}

# Remove prefixes from taxonomy levels
tax_data[, (all_taxa_levels) := lapply(.SD, function(x) str_replace_all(x, "^[a-z]__", "")), .SDcols = all_taxa_levels]

# Remove unknown or NA values from taxonomy
tax_data[, (all_taxa_levels) := lapply(.SD, function(x) str_replace_all(x, "__unknown", "")), .SDcols = all_taxa_levels]

# Prepare taxonomy table
tax_table <- data.frame(tax_data[, ..all_taxa_levels])
rownames(tax_table) <- kraken_data$OTU_ID

# Prepare OTU abundance table
otu_table <- data.frame(kraken_data[, .(OTU_ID, DNA.profile)])
otu_table <- as.data.frame(otu_table[, -1])
rownames(otu_table) <- kraken_data$OTU_ID
colnames(otu_table) <- 'Zymo_D6323'

# Read in the unclassified reads
unclassified_data <- fread("../standards/Zymo.D6323.DNA.report/in2165.DNA/ZymoBIOMICS.fecal.ref.illumina.pe/AbundanceTables/ReadDistributionTable.csv")
setnames(unclassified_data, "#OTU ID", "OTU_ID")

# Add unclassified reads as "unassigned"
unassigned_row <- data.frame(Zymo_D6323 = unclassified_data[OTU_ID == "Unclassified", DNA.profile])
rownames(unassigned_row) <- "unassigned"

# Combine unclassified reads with OTU abundance table
otu_table <- rbind(otu_table, unassigned_row)

tax_table <- rbind(tax_table,
                   data.frame(Kingdom = "unassigned", Phylum = "unassigned", Genus = "unassigned", Species = "unassigned", X = "unassigned", Y = "unassigned", Z = "unassigned", row.names = "unassigned"))

# Create phyloseq object
ps <- phyloseq(otu_table(as.matrix(otu_table), taxa_are_rows = TRUE),
               tax_table(as.matrix(tax_table)))


tax_data <- merge(tax_table, otu_table, by=0)

# Check phyloseq object
print(ps)


fwrite(tax_data, '../standards/Zymo.D6323.DNA.report/zymo_D6323_data.tsv', sep='\t')
fwrite(tax_data[,c("Species", "Zymo_D6323")], '../standards/Zymo.D6323.DNA.report/zymo_D6323_specdata.tsv', sep='\t')

saveRDS(ps, '../standards/Zymo.D6323.DNA.report/zymo_D6323_PS.rds')
