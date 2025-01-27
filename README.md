# Structure
## process_data
- contains notebooks to preprocess datasets and subset for microglia cells, in preparation for integration
- include SEA-AD, SUN-MIT-ROSMAP, Hyman, Rexach, Prater datasets (TODO: add annotation and link to the datasets)
## integration
- scripts to perform integration and save the integrated datasets and models
- `multimil_all.py` integrates MTG and DLPFC parts of SEA-AD and the SUN dataset, excluding the unmatched ATAC-seq part (where there are no matching RNA-seq for the cells).
- `multimil_w_single_atac.py` integrates MTG and DLPFC parts of SEA-AD and the SUN dataset, including all ATAC-seq from SEA-AD dataset.
- `peakvi.py` integrates only the ATAC-seq data from the SEA-AD MTG dataset, was intended to show the unintegratability of matched and unmatched ATAC-seq.
- `scvi_rna.py` integrates only RNA-seq from SEA-AD and SUN dataset
## disease_prediction
- contains notebooks that perform and evaluate patient level predictions based on the integration result.
- `disease_pre.ipynb` shows the results for the integrated dataset obtained from `multimil_all.py`. It includes predictions of disease/healthy, APOE genotype, Thal, and Braak stage. 
- `disease_pre_scvi.ipynb` showa the results for the integrated dataset obtained from `scvi_rna.py`
