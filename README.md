# Heart_MERFISH_analysis

This repository contains the code for the analysis of [Spatially organized cellular communities form the developing human heart](https://www.nature.com/articles/s41586-024-07171-z). 

Raw (count matrices as `.csv` files) and processed data (`.rds` and `.h5ad` files) are accessible at [Dryad](https://datadryad.org/stash/dataset/doi:10.5061/dryad.w0vt4b8vp). 

The notesbooks in [analysis](analysis) folder are used to generate the figures in the manuscript. They should be read in the following order: 
1. [MERFISH_processing](analysis/MERFISH_processing.ipynb)
2. [MERFISH_analysis](analysis/MERFISH_analysis.ipynb)
3. [MERFISH_ventricle_community_analysis](analysis/MERFISH_ventricle_community_analysis.ipynb)
4. [MERFISH_ventricle_depth_analysis](analysis/MERFISH_ventricle_depth_analysis.ipynb)
5. [MERFISH_scRNAseq_integration](analysis/MERFISH_scRNAseq_integration.ipynb)
6. [community_signaling_analysis](analysis/community_signaling_analysis.ipynb)
7. [MERFISH_replicate_heart_analysis](analysis/MERFISH_replicate_heart_analysis.ipynb)