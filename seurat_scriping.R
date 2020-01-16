library(Seurat)
library(dplyr)
read_data <- function(){
  # get list of files with PROJECT/CELL.tsv
  rnaseq_counts_folder <- "/data/gene_count_length_files/"
  rnaseq_counts_projects_names <- list.files(path=rnaseq_counts_folder)
  rnaseq_counts_projects <- vector(mode="list", length=length(rnaseq_counts_projects_names))
  names(rnaseq_counts_projects) <- rnaseq_counts_projects_names
  
  for (project in rnaseq_counts_projects_names){
    full_project_path <- paste(rnaseq_counts_folder, project, sep='')
    project_cell_files <- list.files(path=full_project_path, recursive=TRUE, pattern="*.tsv")
    full_project_df <- data.frame(matrix(ncol = 1, nrow = 0))
    colnames(full_project_df) <- c('Geneid')
    for (file in project_cell_files) {
      # read in each file
      file_df <- read.table(file=paste(full_project_path, file, sep='/'), sep="\t", header=TRUE)
      cell <- strsplit(file, '.tsv')[[1]]
      colnames(file_df) <- gsub(".*\\.bam", cell, colnames(file_df))
      full_project_df <- merge(x=full_project_df, y=file_df, all=TRUE, by='Geneid')
      full_project_df$Length <- coalesce(full_project_df$Length.x, full_project_df$Length.y)
      full_project_df$Length.x <- NULL
      full_project_df$Length.y <- NULL
    }
    rownames(full_project_df) <- full_project_df$Geneid
    rnaseq_counts_projects$project <- full_project_df
  }
}

