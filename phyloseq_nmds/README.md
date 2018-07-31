# Phloseq_NMDS

A Galaxy tool to produce NMDS plots using Phyloseq from either a BIOM1 file or 2 input tables.

Currently produces the plots embedded in a html file for output with links to a PDF file.

Requires:

Phyloseq 1.22.3
r-getopt 1.20.0
ghostscript 9.18


### Run phyloseq_nmds.R with three input files
Rscript phyloseq_nmds.R --otu_table=GP_OTU_TABLE.txt --tax_table=GP_TAX_TABLE.txt --meta_table=GP_SAMPLE_TABLE.txt --method="bray" --kingdom=2 --cutoff=5 --category=6 --outdir=/outputdir --htmlfile=test.html

### Run phyloseq_nmds.R with biom file
Rscript phyloseq_nmds.R --biom=GP.biom --subset=6 --method=NMDS --distance=bray --kingdom=Phylum --cutoff=5 --keep=5 --outdir=/outputdir --htmlfile=biom_out.html

## Version history:

**XML Wrapper:**

Alpha version by Michael Thang of QFAB, Australia.

1.22.3.1: Simon Gladman Melbourne Bioinformatics
    * Incorporated tests
    * Requirements
    * Version statement
    * Citations

1.22.3.2: Michael Thang QFAB, Simon Gladman Melbourne Bioinformatics
    * Uses new version of BIOM1 datatype to get metadata
    * Output label changed as per user requirements


**R Script: phyloseq_nmds.R**

0.1.0: Michael Thang QFAB
    * Original version

0.1.1: Michael Thang QFAB
    * Added extra BIOM import functionality so it doesn't solely rely on phyloseq's internal importer.

0.1.2: Michael Thang QFAB
    * BIOM functionality now requires the column header name in text.
