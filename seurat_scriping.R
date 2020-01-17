library(Seurat)
library(dplyr)

# https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html as a reference

read_data <- function(rnaseq_counts_folder, normalized){
  
  # rnaseq_counts_folder should have a trailing slash
  # normalized indicates if this is jim's normalized data or not
  
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
      if(normalized){
        file_df <- read.table(file=paste(full_project_path, file, sep='/'))
      }else{
        file_df <- read.table(file=paste(full_project_path, file, sep='/'), sep="\t", header=TRUE)
      }
      
      
      # alter column name of hits to be that of the cell
      cell <- strsplit(file, '.tsv')[[1]]
      colnames(file_df) <- gsub(".*\\.bam", cell, colnames(file_df))
      if(normalized){
        file_df$Length <- NULL
        colnames(file_df)[1] <- 'Length'
        file_df$Geneid <- as.numeric(rownames(file_df))
        rownames(file_df) <- NULL
      }
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

write_list_of_dfs <- function(list_of_dfs, dir){
  # dir ends in /
  for(project in names(list_of_dfs)){
    filename <- paste(paste(dir, project, sep=''), ".tsv", sep="")
    print(filename)
    write.table(list_of_dfs[[project]], file=filename)
  }
}

# umbrella function to run pca on all rnaseq counts objects from the projects
run_seurat_pca_on_projects <- function(project_counts, noramlize){
  seurat_objs <- vector(mode="list", length=length(project_counts))
  names(seurat_objs) <- names(project_counts)
  for(project in names(project_counts)){
    pca_dims <- init_seurat_and_run_pca(project_counts[[project]], project, noramlize)
    seurat_objs[[project]] <- pca_dims
  }
  return(seurat_objs)
}

### USE plot_pca FUNCTION TO DETERMINE THE NUMBER OF DIMENSIONS TO PASS INTO CLUSTERING

# run clustering on the entirety
run_cluster_on_seurat_objs <- function(seurat_pca_objs, num_dim_list){
  seurat_cluster_objs <- vector(mode="list", length=length(seurat_pca_objs))
  names(seurat_cluster_objs) <- names(seurat_pca_objs)
  for(project in names(seurat_pca_objs)){
    clusters <- run_cluster(seurat_pca_objs[[project]], num_dim_list[[project]])
    seurat_cluster_objs[[project]] <- clusters
  }
  return(seurat_cluster_objs)
}

# get cluster counts for all
get_cluster_counts <- function(seurat_cluster_objs){
  num_cluser_list <- vector(mode="list", length=length(seurat_cluster_objs))
  names(num_cluser_list) <- names(seurat_cluster_objs)
  for(project in names(seurat_cluster_objs)){
    num_clusters <- get_num_clusters(seurat_cluster_objs[[project]])
    num_cluser_list[[project]] <- num_clusters
  }
  return(num_cluser_list)
}

#get embeddings for all proejcts
get_project_umap_embeddings <- function(seurat_cluster_objs){
  embeddings <- vector(mode="list", length=length(seurat_cluster_objs))
  names(embeddings) <- names(seurat_cluster_objs)
  for(project in names(seurat_cluster_objs)){
    project_embeddings <- get_umap_embeddings(seurat_cluster_objs[[project]])
    embeddings[[project]] <- project_embeddings
  }
  return(embeddings)
  
}

#get assignemnts for all proejcts
get_project_cluster_assignments <- function(seurat_cluster_objs){
  embeddings <- vector(mode="list", length=length(seurat_cluster_objs))
  names(embeddings) <- names(seurat_cluster_objs)
  for(project in names(seurat_cluster_objs)){
    project_embeddings <- get_cluster_assignments(seurat_cluster_objs[[project]])
    embeddings[[project]] <- project_embeddings
  }
  return(embeddings)
  
}

# utility functions for suerat

# 1
init_seurat_and_run_pca <- function (counts_df, project, normalize=TRUE){
  seurat_obj <- CreateSeuratObject(counts_df, proj=project, min.cells = 3, min.features = 200)
  if(normalize){
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  }
  

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

get_num_clusters <- function(seurat_obj){
  
  length(unique(Idents(seurat_obj)))
  
}

# 4
plot_cluster <- function(seurat_obj){
  
  # plot
  
  DimPlot(seurat_obj, reduction = "umap")
  
}

get_umap_embeddings <- function(seurat_obj){
  
  return(seurat_obj@reductions$umap@cell.embeddings)
  
}

get_cluster_assignments <- function(seurat_obj){
  
  return(Idents(seurat_obj))
  
}



# Our Script
'
source("~/single-cell-comparison/seurat_scriping.R")
# rnaseq_counts_projects <- read_data("/data/gene_count_length_files/", FALSE)
load("~/R_Data/counts_data.RData")
pca_seurat_dims <- run_seurat_pca_on_projects(rnaseq_counts_projects, TRUE)
num_dim_list <- vector(mode="list", length=length(pca_seurat_dims))
names(num_dim_list) <- names(pca_seurat_dims)
plot_pca(pca_seurat_dims$SRP011546, TRUE)
num_dim_list$SRP011546 <- 12
plot_pca(pca_seurat_dims$SRP050499, TRUE)
num_dim_list$SRP050499 <- 15
plot_pca(pca_seurat_dims$SRP057196, TRUE)
num_dim_list$SRP057196 <- 12
plot_pca(pca_seurat_dims$SRP061549, TRUE)
num_dim_list$SRP061549 <- 7
plot_pca(pca_seurat_dims$SRP066632, TRUE)
num_dim_list$SRP066632 <- 15
clustering_seurat_objs_counts <- run_cluster_on_seurat_objs(pca_seurat_dims, num_dim_list)
cluster_assignment_list_counts <- get_project_cluster_assignments(clustering_seurat_objs_counts)
cluster_counts_counts <- get_cluster_counts(clustering_seurat_objs_counts)
write_list_of_dfs(cluster_assignment_list_counts, "/data/cluster_assignments_counts/")

load("~/R_Data/normalized_data.RData")


# FALSE for normalization to be performed by Seurat
pca_seurat_dims <- run_seurat_pca_on_projects(normalized_rpkm_data, FALSE)
num_dim_list <- vector(mode="list", length=length(pca_seurat_dims))
names(num_dim_list) <- names(pca_seurat_dims)
plot_pca(pca_seurat_dims$SRP011546, TRUE)
num_dim_list$SRP011546 <- 13
plot_pca(pca_seurat_dims$SRP050499, TRUE)
num_dim_list$SRP050499 <- 12
plot_pca(pca_seurat_dims$SRP057196, TRUE)
num_dim_list$SRP057196 <- 12
plot_pca(pca_seurat_dims$SRP061549, TRUE)
num_dim_list$SRP061549 <- 7
plot_pca(pca_seurat_dims$SRP066632, TRUE)
num_dim_list$SRP066632 <- 17
clustering_seurat_objs_rpkm <- run_cluster_on_seurat_objs(pca_seurat_dims, num_dim_list)
cluster_assignment_list_rpkm <- get_project_cluster_assignments(clustering_seurat_objs_rpkm)
cluster_counts_rpkm <- get_cluster_counts(clustering_seurat_objs_rpkm)
write_list_of_dfs(cluster_assignment_list_rpkm, "/data/cluster_assignments_rpkm/")

pca_seurat_dims <- run_seurat_pca_on_projects(normalized_tpm_data, FALSE)
num_dim_list <- vector(mode="list", length=length(pca_seurat_dims))
names(num_dim_list) <- names(pca_seurat_dims)
plot_pca(pca_seurat_dims$SRP011546, TRUE)
num_dim_list$SRP011546 <- 13
plot_pca(pca_seurat_dims$SRP050499, TRUE)
num_dim_list$SRP050499 <- 12
plot_pca(pca_seurat_dims$SRP057196, TRUE)
num_dim_list$SRP057196 <- 12
plot_pca(pca_seurat_dims$SRP061549, TRUE)
num_dim_list$SRP061549 <- 9
plot_pca(pca_seurat_dims$SRP066632, TRUE)
num_dim_list$SRP066632 <- 17
clustering_seurat_objs_tpm <- run_cluster_on_seurat_objs(pca_seurat_dims, num_dim_list)
cluster_assignment_list_tpm <- get_project_cluster_assignments(clustering_seurat_objs_tpm)
cluster_counts_tpm <- get_cluster_counts(clustering_seurat_objs_tpm)
write_list_of_dfs(cluster_assignment_list_tpm, "/data/cluster_assignments_tpm/")
'