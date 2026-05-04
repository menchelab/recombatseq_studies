### Generates a simulated experimental design matrix with confounded variables
createConfoundedDesign <- function(n_samples) {
  new_group <- replicate(3, sample(0:2, n_samples, replace=TRUE))
  last_col <- (new_group[,1]+1)%%3
  covmat <- data.frame()
  covmat <- cbind(new_group, last_col)
  covmat <- as.data.frame(covmat)
  covmat[] <- lapply(covmat, as.factor)
  
  return(covmat)
}

### edgeR-based differential expression pipeline, taken from ComBat-seq
edgeR_DEpipe <- function(counts_mat, batch, group, include.batch, alpha.unadj, alpha.fdr, covar=NULL){
  cat("DE tool: edgeR\n")
  y <- DGEList(counts=counts_mat)
  y <- edgeR::calcNormFactors(y, method="TMM")
  if(include.batch){
    cat("Including batch as covariate\n")
    design <- model.matrix(~ as.factor(group) + as.factor(batch))
  }else{
    cat("Default group as model matrix\n")
    design <- model.matrix(~as.factor(group))
  }
  if(!is.null(covar)){
    cat("Including surrogate variables or unwanted variation variables\n")
    design <- cbind(design, covar)
  }
  y <- edgeR::estimateDisp(y, design)
  fit <- edgeR::glmQLFit(y, design)
  qlf <- edgeR::glmQLFTest(fit, coef=2)
  de_res <- edgeR::topTags(qlf, n=nrow(counts_mat))$table
  
  de_called <- rownames(de_res)[de_res$PValue < alpha.unadj]
  de_called_fdr <- rownames(de_res)[de_res$FDR < alpha.fdr]
  return(list(unadj=de_called, fdr=de_called_fdr, de_res=de_res, design=design))
}

### Computes performance metrics (TPR, FPR, precision) for detected DE genes against a ground truth set
perfStats <- function(called_vec, ground_truth_vec, N_genes){
  if(length(called_vec)==0){
    tpr <- fpr <- 0
    prec <- NA
  }else{
    tp <- length(intersect(called_vec, ground_truth_vec))
    fp <- length(setdiff(called_vec, ground_truth_vec))
    N_DE <- length(ground_truth_vec)
    N_nonDE <- N_genes - N_DE
    
    tpr <- tp / N_DE
    fpr <- fp / N_nonDE
    prec <- tp / length(called_vec)
  }
  return(c(tpr=tpr, fpr=fpr, prec=prec))
}
