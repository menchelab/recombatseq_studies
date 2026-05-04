library(airway)
library(reComBatseq)
library(splatter)
library(stringr)
source("../helper_functions.R")

# Initialize base parameters
data(airway)
airway_counts <- as.matrix(assay(airway))
params <- splatEstimate(airway_counts)
params <- newSplatParams()

# Set up experiment parameters
folder <- "batch_study"
n_batches <- 24

# Fixed params
gene_count <- 2^12
sample_count <- 6000   # 2000 * n_batches
defacLoc = 0.5
defacScale = 0.2

# Create directory if not available
path_base <- paste0(folder,"/experiment_",as.character(n_batches))
dir.create(path_base, showWarnings = FALSE)

# Simulation pipeline
res_list <- vector("list", 5)

for(iter in 1:5){
  cat(paste("Iteration", iter, "\n"))
  seed <- sample.int(1e6, 1)
  facLoc = abs(rnorm(n_batches, mean=2, sd=0.3))  
  facScale = abs(rnorm(n_batches, mean=0.8, sd=0.15)) 
  
  
  sim_batch <- splatSimulate(nGenes=gene_count, batchCells=rep(sample_count/n_batches, n_batches),
                             mean.rate = params@mean.rate, mean.shape = params@mean.shape,
                             bcv.common = params@bcv.common, bcv.df = params@bcv.df,
                             group.prob=c(0.5,0.5), de.prob=0.1,
                             batch.facLoc = facLoc, batch.facScale = facScale,
                             de.facLoc = defacLoc, de.facScale = defacScale,
                             method='groups',
                             verbose=F,
                             batch.rmEffect=FALSE, seed=seed)
  
  sim_noBatch <- splatSimulate(nGenes=gene_count, batchCells=rep(sample_count/n_batches,n_batches),
                               mean.rate = params@mean.rate, mean.shape = params@mean.shape,
                               bcv.common = params@bcv.common, bcv.df = params@bcv.df,
                               group.prob=c(0.5,0.5), de.prob=0.1,
                               batch.facLoc = facLoc, batch.facScale = facScale,
                               de.facLoc = defacLoc, de.facScale = defacScale,
                               method='groups',
                               verbose=F,
                               batch.rmEffect=TRUE, seed=seed)
  
  batch <- colData(sim_batch)@listData[["Batch"]]
  group <- colData(sim_batch)@listData[["Group"]]
  batch <- as.factor(as.numeric(unlist(str_extract_all(batch, "\\d+"))))
  group <- as.factor(as.numeric(unlist(str_extract_all(group, "\\d+"))))
  
  metadata <- cbind(batch, group)
  write.csv(metadata,
            paste0(path_base,"/iter",iter,"_metadata_BS.csv"))
  
  batch_df <- as.matrix(counts(sim_batch))
  nobatch_df <- as.matrix(counts(sim_noBatch))
  
  write.csv(batch_df,
            paste0(path_base,"/iter",iter,"_batch_df_BS.csv"))
  write.csv(nobatch_df,
            paste0(path_base,"/iter",iter,"_nobatch_df_BS.csv"))
  
  count_batch_transformed <- DGEList(counts=batch_df)
  count_batch_transformed <- edgeR::calcNormFactors(count_batch_transformed, method="TMM")
  count_batch_transformed <- voom(count_batch_transformed, model.matrix(~as.factor(group)))
  
  write.csv(count_batch_transformed$E,
            paste0(path_base,"/iter",iter,"_countmat_batch_transformed_BS.csv"))
  
  # Create confounded design matrix
  covmat <- createConfoundedDesign(sample_count)
  covmat$group <- group
  #qr(covmat)$rank
  
  write.csv(covmat,
            paste0(path_base,"/covmat_BS.csv"))
  
  # Batch correction
  start.time <- Sys.time()
  recombatseq_df <- reComBat.seq(batch_df, batch = batch, wanted.variation = covmat, num.threads=2)
  end.time <- Sys.time()
  time_reg <- as.numeric(difftime(end.time,start.time, units="mins"))
  
  write.csv(recombatseq_df,
            paste0(path_base,"/iter",iter,"_recombatseq_df_BS.csv"))

  # DE Analysis
  de_cols <- grep("^DEFacGroup", colnames(rowData(sim_batch)), value = TRUE)
  de_sums <- rowSums(as.matrix(rowData(sim_batch)[, de_cols]))
  de_ground_truth <- rownames(rowData(sim_batch))[de_sums != length(de_cols)]
  
  saveRDS(de_ground_truth, 
          file = paste0(path_base,"/iter",iter,"_DEgenes_BS.rds"))
  
  nobatch_cor <- edgeR_DEpipe(nobatch_df, batch=batch, group=group,
                              include.batch=FALSE, alpha.unadj=0.05, alpha.fdr=0.01)[["fdr"]]
  batch_cor <- edgeR_DEpipe(batch_df, batch=batch, group=group,
                            include.batch=FALSE, alpha.unadj=0.05, alpha.fdr=0.01)[["fdr"]]
  recombatseq_cor <- edgeR_DEpipe(recombatseq_df, batch=batch, group=group,
                                   include.batch=FALSE, alpha.unadj=0.05, alpha.fdr=0.01)[["fdr"]]
  
  stats_nobatch <- perfStats(nobatch_cor, de_ground_truth, gene_count)
  stats_batch <- perfStats(batch_cor, de_ground_truth, gene_count)
  stats_recombatseq <- perfStats(recombatseq_cor, de_ground_truth, gene_count)
  
  # RESULTS DATA FRAME
  res_df <- data.frame(
    time=c(0, 0, time_reg, 0),
    tpr=c(stats_batch[["tpr"]], stats_nobatch[["tpr"]], stats_recombatseq[["tpr"]], 0),
    fpr=c(stats_batch[["fpr"]], stats_nobatch[["fpr"]], stats_recombatseq[["fpr"]], 0),
    prec=c(stats_batch[["prec"]], stats_nobatch[["prec"]], stats_recombatseq[["prec"]], 0)
  )
  
  res_df <- t(res_df)
  colnames(res_df) <- c("withBatch", "withoutBatch", "reComBat-seq", "reComBat")
  
  res_list[[iter]] <- res_df
}

res_df <- Reduce(`+`, res_list) / length(res_list)
res_df

#### save results
write.csv(res_df,
          paste0(path_base,"/results_BS.csv"))

