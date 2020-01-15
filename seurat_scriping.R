library(Seurat)

# get list of files with PROJECT/CELL.tsv
rnaseq_counts_folder <- "/data/gene_count_length_files/"
rna_seq_counts_files <- list.files(path=rnaseq_counts_folder, recursive=TRUE, pattern="*.tsv")
full_df <- data.frame()

for (file in rna_seq_counts_files) {
  # read in each file
  file_df <- read.table(file=paste(rnaseq_counts_folder, file, sep=''), sep="\t", header=TRUE, row.names='Geneid')
  file_name <- strsplit(file, '/')[[1]][2]
  cell <- strsplit(file_name, '.tsv')[[1]]
  colnames(file_df) <- gsub(".*\\.bam", cell, colnames(file_df))
  full_df <- merge(x=full_df, y=file_df, all=TRUE)
  
}