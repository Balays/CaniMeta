# CaniMeta

Analysis and quality-control workflow for our canine fecal microbiome dataset, combining Illumina and Oxford Nanopore WGS, extraction-method comparisons and full-length 16S primer testing.

**Paper:** [Canine Fecal Microbiome Dataset: Ultra-deep Multi-platform Sequencing Across Extraction and Library Protocols](https://doi.org/10.1038/s41597-026-07594-5)

The main workflow is `CaniMeta.Rmd`. It processes the ENA datasets, summarizes read length and quality, compares sequencing and extraction protocols, and generates the figures and tables used in the paper.

## Repository structure

- `1_Dog.M0_WGS/` – ultra-deep Illumina and ONT WGS comparison
- `2_40dogs_ZHMW_vs_ZMB/` – longitudinal comparison of two DNA extraction methods
- `3_Primer_comparison/` – full-length 16S primer comparison
- `scripts/` – download, quality-control and helper scripts
- `ENA/` – accession metadata and read-level QC summaries

Some preprocessing chunks in `CaniMeta.Rmd` are set to `eval=FALSE` because they only need to be run once. The paths to local helper functions near the top of the file should be changed before running the complete workflow.
