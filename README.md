# Using scSRASeq to verify and improve submitted metadata availble in SRA

## Team members
 - [Jim Abraham](https://github.com/jcabraham)
 - [Ilan Gold](https://github.com/ilan-gold)
 - [Vamsi Kodali](https://github.com/vkkodali)
 - [Lukas Wagner](https://github.com/lwagnerdc)

## The Problem
Identifying publicly submitted projects relevant to a research goal for baselines or benchmarking is often desirable.
 However, the usability and accuracy of metadata for RNA-Seq data is variable.  Methods for validating available metadata or potentially characterizing datasets where clear characterization in metadata is lacking would be a useful resource and also a valuable exercise for codeathon participants.
 Internal consistency, comparison with selected reference sets (either bulk or other single-cell), and aggregate properties of alignment against genome can all be used for this comparison. Alignment against genome could be based on methods in use now, as could comparison with selected reference datasets for which some results have been precomputed.

## Methods and Approach
We will use basic data processing tools like `Python`, shell-scripting, and `R` to do data-wrangling. We will use a variety of statistical methods (some basic, like [Spearman Rank Correlation](https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient), and some bioinformatics-specific, like [Seurat](https://satijalab.org/seurat/)) in `R` to assess accuracy of metadata submitted for fields like cell type and number of cells. For example, we can compare the output of Seurat's clustering (i.e the number of clusters) against what is submitted in the metadata.  A large deviation would suggest something is wrong with the metadata.

More specifically, we downloaded the metadata, tabulated it using `Perl` and shell scripting (see [Tabluating_Metadata.md](./Tabulating_Metadata.md)).  This gives us a rough way of looking at "cell types."  We then took cleaned gene expression data aligned with Hisat2 and passed it through a Seurat pipeline similar to the one on the [website](https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html).
Hisat2 was configured to retain up to 5 alignments against GRCh38.  FeatureCounts was used with these alignments and NCBI's annotation of the genome to generate gene counts.  Options to include multiply placed reads and overlapping features were used, with fractional weighting for such reads in the gene feature counts.  The output of Seurat, the clusters, should represent cell types and was compared to the cell types obtained from downloading the metadata.

## Impact
Using metadata and generating metrics for comparison can help to identify interesting features in research results.
Comparing public data with results generated in the course of research can improve reproducibility and robustness of computational tools used for research analysis. Having researchers identify possible improvements to existing metadata could improve submission practice as well as disseminate knowledge of best practices for using existing metadata in SRA.

## Datasets Used
We downloaded [SRA](https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/) files.  SRA single cell projects which met the following criteria were used as illustrative examples:
- mostly human
- one to one run to cell mapping for most of the project
- 100 to 10000 cells sampled
- interesting range of sample attributes (disease or developmental ES stages)

## [Results](./RESULTS.md)
