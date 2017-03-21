# EMPHASIS EPIC analysis
## 1. Methylation data processing
1. Run `asaffa/EMPHASIS/EPIC_analysis/01_meffil_QC.R path cores` as an R batch job where:
* `path` = path to the project folder which should contain the .idat files as well as the sample sheet (.csv)
* `cores` = number of cores to use for processing, default is 8
2. Inspect output from the script to ensure it completed running succesfully
3. Examine the generated QC reports (in `path/results/`)
4. Look at `path/pc_fit.pdf` and determine elbow point of the plot, this will be the number of PCs used for functional normalization in the next step
5. Run `asaffa/EMPHASIS/EPIC_analysis/02_meffil_norm.R path PCs batch cores` as an R batch job where:
* `path` = as above
* `PCs` = number of PCs to use (see 4.)
* `batch` = list of batch variables to adjust for, this should be a single row .tsv file in the project dir
* `cores` = as above
6. Inspect output from the script
7. Examine the generated normalization reports (in `path/results/`)

Example code:
* Rscript 01_meffil_QC.R /home/EMPHASIS 16
* Rscript 02_meffil_norm.R /home/EMPHASIS 14 batch.tsv 16

## 2. Outcome data processing
## 3. Intervention-methylation associations
## 4. Methylation-outcome associations
## 5. Pathways analysis
## 6. Cross-cohort comparison
## 7. Causal inference and mQTL analysis
