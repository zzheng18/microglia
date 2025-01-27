import numpy as np
import scanpy as sc
import scvi
import torch

adata = sc.read_h5ad('/home/icb/zihe.zheng/projects/microglia/data/SEAAD/microglia_SEAAD_MTG_ATAC.h5ad')

adata = adata[~adata.obs['Supertype'].isin(['n_genes', 'fraction_mito', 'doublet_score','cluster_fraction_mito_flag','cluster_doublet_score_flag','cluster_fraction_ribo_flag', 'cluster_n_genes_flag'])]

print("# regions before filtering:", adata.shape[-1])

# compute the threshold: 5% of the cells
min_cells = int(adata.shape[0] * 0.05)
# in-place filtering of regions
sc.pp.filter_genes(adata, min_cells=min_cells)

print("# regions after filtering:", adata.shape[-1])

scvi.model.PEAKVI.setup_anndata(adata)

model = scvi.model.PEAKVI(adata)
model.train()

model.save('/home/icb/zihe.zheng/projects/microglia/peakvi', overwrite=True)

PEAKVI_LATENT_KEY = "X_peakvi"

latent = model.get_latent_representation()
adata.obsm[PEAKVI_LATENT_KEY] = latent
print(latent.shape)

adata.obs = adata.obs.astype(str)
adata.write_h5ad('/home/icb/zihe.zheng/projects/microglia/data/integrated_atac.h5ad')

