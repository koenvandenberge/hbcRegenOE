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