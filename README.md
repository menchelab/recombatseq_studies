# reComBat-seq Studies
Summary of the pipelines and experiments associated with `reComBat-seq`

## File descriptions

### Simulation Studies

#### Gene Count Study

+ **gene_study_step1.R**
    + Simulates RNA-seq datasets with and without batch effects - 4,096 samples and genes varying between 128 and 16,384.
    + Applies `ComBat-seq` and `reComBat-seq` for batch correction, runs DE analysis (edgeR) on all data, computes precision/TPR/FPR.
    + Saves raw counts, normalized counts (voom/TMM), confounded covariate matrices and corrected count matrices as CSV files.
+ **gene_study_step2.ipynb**
    + Applies `reComBat` for batch correction on the simulated data from the previous step, using the transformed normalized counts instead of the raw matrix.
+ **gene_study_step3.R**
    + Runs DE analysis on the `reComBat` (limma) corrected data.

#### Sample Count Study

+ **sample_study_step1.R**
    + Simulates RNA-seq datasets with and without batch effects - 4,096 genes and samples varying between 16 and 16 384.
    + Applies `ComBat-seq` and `reComBat-seq` for batch correction, runs DE analysis (edgeR) on all data, computes precision/TPR/FPR.
    + Saves raw counts, normalized counts (voom/TMM), confounded covariate matrices and corrected count matrices as CSV files.
+ **sample_study_step2.ipynb**
    + Applies `reComBat` for batch correction on the simulated data from the previous step, using the transformed normalized counts instead of the raw matrix.
+ **results_pipeline_step3.R**
    + Runs DE analysis on the `reComBat` (limma) corrected data.
 
#### Batch Count Study

+ **batch_study_step1.R**
    + Simulates RNA-seq datasets with and without batch effects for two design scenarios - 6000 samples/ 2000 samples per batch, 4,096 genes and batches varying between 2 and 40.
    + Applies `reComBat-seq` for batch correction, runs DE analysis (edgeR) on all data, computes precision/TPR/FPR.
    + Saves raw counts, normalized counts (voom/TMM), confounded covariate matrices and corrected count matrices as CSV files.
+ **batch_study_step2.ipynb**
    + Applies `reComBat` for batch correction on the simulated data from the previous step, using the transformed normalized counts instead of the raw matrix.
+ **batch_study_step3.R**
    + Runs DE analysis on the `reComBat` (limma) corrected data.

 
### Real Data
Code to reproduce results regarding algorithm comparisons on real data and other results shown in the accompanying paper can be found in `algo_comparisons`. Data is available on [Zenodo](https://doi.org/10.5281/zenodo.19736515)
