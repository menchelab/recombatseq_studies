library(limma)
source("../helper_functions.R")

# Set up experiment parameters
folder <- "gene_study"
gene_count <- 2^7
experiment <- "2(7)"

# Fixed params
path_base <- paste0(folder,"/experiment_",experiment)
n_batches <- 2
sample_count <- 2^12

stats_recombat <- vector("list", 5)

for(iter in 1:5){
  cat(paste("Iteration", iter, "\n"))
  
  recombat_df <- read.csv(paste0(path_base, "_iter", iter,"_recombat_df_GS.csv"), row.names=1)
  metadata <- read.csv(paste0(path_base, "_iter", iter,"_metadata_GS.csv"), row.names=1)
  group <- as.factor(metadata[["group"]])
  
  
  # Differential expression
  de_genes <- readRDS(paste0(path_base, "_DEgenes_GS.rds"))
  de_genes <- as.numeric(unlist(str_extract_all(de_genes, "\\d+")))
  
  # reComBat
  vfit <- lmFit(scale(recombat_df), model.matrix(~as.factor(group)))
  efit <- eBayes(vfit)
  tests <- decideTests(efit)
  recombat_cor <- which(tests[,2]!=0)
  
  # Stats
  stats_recombat[[iter]] <- perfStats(recombat_cor, de_genes, gene_count)
}

colMeans(do.call(rbind, stats_recombat))