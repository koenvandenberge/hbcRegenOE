# hbcRegenOE

## Analysis scripts

Below is a brief summary of what each of the scripts do.

 - `01_20200415_preprocessingTI`: Preprocessing of scRNA-seq data to prepare it for trajectory inference. Includes removal of doublets, temporally asynchronous cells as well as respiratory and microvillous cells.
 - `02_20200819_trajectoryInference`: Trajectory inference.
 - `02Supp_geneExpressionTrajectory`: Supplementary script visualizing gene expression of genes along the trajectory or pseudotime. Use this if you'd like to explore gene expression, or the Shiny app at https://scf.berkeley.edu/shiny/koenvdberge/3dOES/.
 - `03_evaluateKHBCRegen`: Evaluating the number of knots to use for `tradeSeq`.
 - `04_fitNBGam_differentKnots`: Fit `tradeSeq` NB-GAM.
 - `05_20201001_figure1`: Plots of trajectory and marker discovery for each lineage using `tradeSeq`. Includes code for generating panels b-e of Figure 1.
 - `06_20200916_hbcInjuryTFAnalysis_threshold`: Transcription factor peak analysis finding TFs that peak in each lineage. Code for Figure 2.
 - `07_network2Core`: Extracting the 2-core of the association networks estimated using *learn2count*.
 - `08_20201007_tfCascadeNetworkIntegration_figure2`: Relate association network and TF peak analyses. Generation of Figure 3.
 - `08Supp_nichenetAnalysis`: Supplementary script on NicheNet analysis.
 - `09_20200610_scATAC_injury_archr_samples_pairwiseHBC_v3`: scATAC-seq data analysis. Identification of cell types, and HBC cell states. Peak calling on HBC scATAC-seq cells. TF motif enrichment in marker peaks for each HBC cell state. Plots for Figure 5.
 - `10_20200715_downstreamArchRAnalysis_v2`: Differential accessibility analysis of `ArchR` gene activity scores. Comparing injured vs uninjured cells over all cell states. Assessing injury effect within hybrid cells. Assessing cell state effect within injured and uninjured cells.
 - `10Supp_atacHBCMotifs_homer`: Homer TF enrichment of DA genes.
 - `11_20201012_ATAC_hbcClustering_v3`: Integration of 24H PI scRNA-seq HBC cells with scATAC-seq HBC cells. Cell type transfer from scRNA-seq to scATAC-seq. Checking gene expression of scATAC-seq markers as well as contribution of injured versus uninjured cells. Identification of scATAC-seq markers for transferred hybrid1 and hybrid2 labels. Peak calling using transferred labels and TF motif enrichment on marker peaks identified as such. Plots for Figure 6.
 - `11Supp_ATACExploration`: Identification of marker genes for each of the hybrid transferred labels. Plotting accessibility of key lineage-specific genes from other cell types and lineages. Plotting accessibility of hybrid markers that were identified using scRNA-seq.
 - `11Supp2_atacRNAExploration`: Exploring ATAC markers in scRNA-seq data. This includes the current lineage-traced data and also (a subset of) the whole OE data.
 - `11Supp3_integratescRNAseq`: Integration of lineage-traced with whole OE data. Shows a similarity between regenerated HBCs from lineage-traced data to truly resting cells from whole OE data. Identification of marker genes of hybrid clusters and activated/resting HBCs using scRNA-seq data.
 - `12_211022_integrateCascadeMotifs`: Supplementary script on linking the results of TF cascade to motif enrichment.

## Supplementary Data


## Software versions

```
## R version 3.6.2 (2019-12-12)
## Platform: x86_64-apple-darwin15.6.0 (64-bit)
## Running under: macOS Mojave 10.14.5
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
##  [1] grid      parallel  stats4    stats     graphics  grDevices utils    
##  [8] datasets  methods   base     
## 
## other attached packages:
##  [1] chromVARmotifs_0.2.0               motifmatchr_1.8.0                 
##  [3] circlize_0.4.8                     ComplexHeatmap_2.3.1              
##  [5] presto_1.0.0                       nabor_0.5.0                       
##  [7] BSgenome.Mmusculus.UCSC.mm10_1.4.0 BSgenome_1.54.0                   
##  [9] rtracklayer_1.46.0                 Biostrings_2.54.0                 
## [11] XVector_0.26.0                     ggalluvial_0.12.0                 
## [13] ArchR_0.9.5                        magrittr_1.5                      
## [15] rhdf5_2.30.1                       data.table_1.12.8                 
## [17] patchwork_1.0.0                    Seurat_3.1.5                      
## [19] Signac_0.2.5                       harmony_1.0                       
## [21] Rcpp_1.0.6                         scater_1.14.6                     
## [23] uwot_0.1.5                         scran_1.14.5                      
## [25] forcats_0.4.0                      stringr_1.4.0                     
## [27] dplyr_1.0.3                        purrr_0.3.3                       
## [29] readr_1.3.1                        tidyr_1.0.2                       
## [31] tibble_2.1.3                       tidyverse_1.3.0                   
## [33] irlba_2.3.3                        Matrix_1.3-2                      
## [35] gridExtra_2.3                      wesanderson_0.3.6                 
## [37] pheatmap_1.0.12                    ggplot2_3.3.2                     
## [39] RColorBrewer_1.1-2                 clusterExperiment_2.6.1           
## [41] bigmemory_4.5.36                   rgl_0.100.30                      
## [43] cowplot_1.0.0                      SingleCellExperiment_1.8.0        
## [45] SummarizedExperiment_1.16.1        DelayedArray_0.12.2               
## [47] BiocParallel_1.20.1                matrixStats_0.56.0                
## [49] Biobase_2.46.0                     GenomicRanges_1.38.0              
## [51] GenomeInfoDb_1.22.0                IRanges_2.20.2                    
## [53] S4Vectors_0.24.3                   BiocGenerics_0.32.0               
## [55] slingshot_1.4.0                    princurve_2.1.4                   
## 
## loaded via a namespace (and not attached):
##   [1] R.methodsS3_1.7.1           pkgmaker_0.31              
##   [3] bit64_0.9-7                 knitr_1.28                 
##   [5] R.utils_2.9.2               KEGGREST_1.26.1            
##   [7] TFBSTools_1.22.0            RCurl_1.98-1.1             
##   [9] doParallel_1.0.15           generics_0.0.2             
##  [11] RSQLite_2.2.0               VGAM_1.1-2                 
##  [13] RANN_2.6.1                  future_1.16.0              
##  [15] ggfittext_0.9.0             bit_1.1-15.2               
##  [17] phylobase_0.8.6             webshot_0.5.2              
##  [19] xml2_1.3.2                  lubridate_1.7.4            
##  [21] httpuv_1.5.2                assertthat_0.2.1           
##  [23] DirichletMultinomial_1.26.0 viridis_0.5.1              
##  [25] xfun_0.12                   hms_0.5.3                  
##  [27] evaluate_0.14               promises_1.1.0             
##  [29] fansi_0.4.1                 progress_1.2.2             
##  [31] caTools_1.18.0              dbplyr_1.4.2               
##  [33] readxl_1.3.1                igraph_1.2.4.2             
##  [35] DBI_1.1.0                   htmlwidgets_1.5.1          
##  [37] ellipsis_0.3.0              gggenes_0.4.0              
##  [39] RSpectra_0.16-0             crosstalk_1.0.0            
##  [41] backports_1.1.5             annotate_1.64.0            
##  [43] gridBase_0.4-7              locfdr_1.1-8               
##  [45] RcppParallel_4.4.4          vctrs_0.3.6                
##  [47] Cairo_1.5-10                here_0.1                   
##  [49] ROCR_1.0-7                  withr_2.1.2                
##  [51] sctransform_0.3.2           GenomicAlignments_1.22.1   
##  [53] prettyunits_1.1.1           softImpute_1.4             
##  [55] cluster_2.1.0               seqLogo_1.50.0             
##  [57] ape_5.3                     lazyeval_0.2.2             
##  [59] crayon_1.3.4                genefilter_1.68.0          
##  [61] labeling_0.3                edgeR_3.28.0               
##  [63] pkgconfig_2.0.3             nlme_3.1-142               
##  [65] vipor_0.4.5                 rlang_0.4.10               
##  [67] globals_0.12.5              lifecycle_0.2.0            
##  [69] miniUI_0.1.1.1              registry_0.5-1             
##  [71] bigmemory.sri_0.1.3         modelr_0.1.5               
##  [73] rsvd_1.0.2                  cellranger_1.1.0           
##  [75] rprojroot_1.3-2             lmtest_0.9-37              
##  [77] rngtools_1.5                ggseqlogo_0.1              
##  [79] Rhdf5lib_1.8.0              zoo_1.8-7                  
##  [81] reprex_0.3.0                beeswarm_0.2.3             
##  [83] GlobalOptions_0.1.1         ggridges_0.5.2             
##  [85] rjson_0.2.20                png_0.1-7                  
##  [87] viridisLite_0.3.0           bitops_1.0-6               
##  [89] R.oo_1.23.0                 rncl_0.8.3                 
##  [91] KernSmooth_2.23-16          blob_1.2.1                 
##  [93] DelayedMatrixStats_1.8.0    shape_1.4.4                
##  [95] manipulateWidget_0.10.0     zinbwave_1.11.6            
##  [97] CNEr_1.20.0                 scales_1.1.0               
##  [99] memoise_1.1.0               plyr_1.8.5                 
## [101] ica_1.0-2                   howmany_0.3-1              
## [103] gplots_3.0.1.2              bibtex_0.4.2.2             
## [105] gdata_2.18.0                zlibbioc_1.32.0            
## [107] compiler_3.6.2              lsei_1.2-0                 
## [109] dqrng_0.2.1                 clue_0.3-57                
## [111] fitdistrplus_1.0-14         Rsamtools_2.2.1            
## [113] cli_2.0.1                   ade4_1.7-13                
## [115] listenv_0.8.0               pbapply_1.4-2              
## [117] MASS_7.3-51.4               tidyselect_1.1.0           
## [119] stringi_1.4.5               yaml_2.2.1                 
## [121] BiocSingular_1.2.1          locfit_1.5-9.1             
## [123] ggrepel_0.8.1               tools_3.6.2                
## [125] future.apply_1.4.0          rstudioapi_0.11            
## [127] uuid_0.1-2                  TFMPvalue_0.0.8            
## [129] foreach_1.4.7               RNeXML_2.4.2               
## [131] farver_2.0.3                Rtsne_0.15                 
## [133] digest_0.6.27               FNN_1.1.3                  
## [135] shiny_1.4.0                 broom_0.5.4                
## [137] later_1.0.0                 RcppAnnoy_0.0.16           
## [139] httr_1.4.1                  AnnotationDbi_1.48.0       
## [141] npsurv_0.4-0                kernlab_0.9-29             
## [143] colorspace_1.4-1            rvest_0.3.5                
## [145] XML_3.99-0.3                fs_1.3.1                   
## [147] reticulate_1.14             splines_3.6.2              
## [149] statmod_1.4.33              plotly_4.9.1               
## [151] xtable_1.8-4                poweRlaw_0.70.2            
## [153] jsonlite_1.6.1              R6_2.4.1                   
## [155] pillar_1.4.3                htmltools_0.5.1.1          
## [157] mime_0.9                    NMF_0.21.0                 
## [159] glue_1.4.2                  fastmap_1.0.1              
## [161] BiocNeighbors_1.4.1         codetools_0.2-16           
## [163] tsne_0.1-3                  lattice_0.20-38            
## [165] ggbeeswarm_0.6.0            leiden_0.3.4               
## [167] gtools_3.8.1                GO.db_3.10.0               
## [169] survival_3.1-8              limma_3.42.1               
## [171] rmarkdown_2.1               munsell_0.5.0              
## [173] GetoptLong_0.1.8            GenomeInfoDbData_1.2.2     
## [175] iterators_1.0.12            HDF5Array_1.14.1           
## [177] haven_2.2.0                 reshape2_1.4.3             
## [179] gtable_0.3.0
```
