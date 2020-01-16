# Using scSRASeq to verify and improve submitted metadata availble in SRA

## Team members
 - Vamsi
 - [Ilan Gold](https://github.com/ilan-gold)
 - Jim Abraham
 - Lukas Wagner

## The Problem
Identifying publicly submitted projects relevant to a research goal for baselines or benchmarking is often desirable. 
 However, the usability and accuracy of metadata for RNA-Seq data is variable.  Methods for validating available metadata or potentially characterizing datasets where clear characterization in metadata is lacking would be a useful resource and also a valuable exercise for codeathon participants. 
 Internal consistency, comparison with selected reference sets (either bulk or other single-cell), and aggregate properties of alignment against genome can all be used for this comparison. Alignment against genome could be based on methods in use now, as could comparison with selected reference datasets for which some results have been precomputed. 
 
 ## Methods and Proposed Solution
 We will use basic data processing tools like Python, shell-scripting, and R to do data-wrangling. We will use a variety of statistical methods (some basic, like [Spearman Rank Correlation](https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient), and some bioinformatics-specific, like [Seurat](https://satijalab.org/seurat/)) in R to assess accuracy of metadata submitted for fields like cell type and number of cells.
 
 ## Impact 
 Using metadata and generating metrics for comparison can help to identify interesting features in research results.
 Comparing public data with results generated in the course of research can improve reproducibility and robustness of computational tools used for research analysis. Having researchers identify possible improvements to existing metadata could improve submission practice as well as disseminate knowledge of best practices for using existing metadata in SRA.
 
 ## Lessons Learned

## Datasets Used
