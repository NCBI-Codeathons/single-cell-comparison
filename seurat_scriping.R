library(Seurat)
library(dplyr)
read_data <- function(rnaseq_counts_folder){
  
  # get list of files with PROJECT/CELL.tsv
  rnaseq_counts_projects_names <- list.files(path=rnaseq_counts_folder)
  rnaseq_counts_projects <- vector(mode="list", length=length(rnaseq_counts_projects_names))
  names(rnaseq_counts_projects) <- rnaseq_counts_projects_names
  
  for (project in rnaseq_counts_projects_names){
    
    # get all cell files in project folder and create empty large df
    full_project_path <- paste(rnaseq_counts_folder, project, sep='')
    project_cell_files <- list.files(path=full_project_path, recursive=TRUE, pattern="*.tsv")
    full_project_df <- data.frame(matrix(ncol = 1, nrow = 0))
    colnames(full_project_df) <- c('Geneid')
    
    for (file in project_cell_files) {
      
      # read in each file
      file_df <- read.table(file=paste(full_project_path, file, sep='/'), sep="\t", header=TRUE)
      
      # alter column name of hits to be that of the cell
      cell <- strsplit(file, '.tsv')[[1]]
      colnames(file_df) <- gsub(".*\\.bam", cell, colnames(file_df))
      
      # merge into full dataframe for the project
      full_project_df <- merge(x=full_project_df, y=file_df, all=TRUE, by='Geneid')
      full_project_df$Length <- coalesce(full_project_df$Length.x, full_project_df$Length.y)
      full_project_df$Length.x <- NULL
      full_project_df$Length.y <- NULL
      
    }
    print(full_project_df)
    #put the project into the project list
    rownames(full_project_df) <- full_project_df$Geneid
    full_project_df$Geneid <- NULL
    rnaseq_counts_projects[[project]] <- full_project_df
    
  }
  return(rnaseq_counts_projects)
}

# 1
init_seurat_and_run_pca <- function (counts_df, project){
  seurat_obj <- CreateSeuratObject(counts_df, proj=project, min.cells = 3, min.features = 200)
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

  # find features
  
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
 
  # scale data
  
  all.genes <- rownames(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj, features = all.genes)
  
  # run PCA
  
  seurat_obj <- RunPCA(seurat_obj)

}
# 2
plot_pca <- function(seurat_obj, quick_run){
  if(quick_run){
    ElbowPlot(seurat_obj)
  } else {
    seurat_obj <- JackStraw(seurat_obj, num.replicate = 100)
    seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:20)
    JackStrawPlot(seurat_obj, dims = 1:15)
  }  
}
# 3
run_cluster <- function(seurat_obj, num_dim) {
    
    # get clusters
  
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:num_dim)
    seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:num_dim)
    
    return(seurat_obj)
}
# 4
plot_cluster <- function(seurat_obj){
  
  # plot
  
  DimPlot(seurat_obj, reduction = "umap")
}