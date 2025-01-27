import numpy as np
import scanpy as sc
import anndata as ad
import scvi
import torch

rna_mtg = sc.read_h5ad('/home/icb/zihe.zheng/projects/microglia/data/SEAAD/reduced_microglia_SEAAD_MTG_RNAseq.h5ad')
rna_dlpfc = sc.read_h5ad('/home/icb/zihe.zheng/projects/microglia/data/SEAAD/reduced_microglia_SEAAD_DLPFC_RNA.h5ad')
rna_sun = sc.read_h5ad('/home/icb/zihe.zheng/projects/microglia/data/SUN/reduced_microglia_MIT_ROSMAP_RNA.h5ad')

rna_sun.layers['UMIs'] = rna_sun.X.copy()

rna_mtg.obs['batch'] = 'seaad_mtg'
rna_dlpfc.obs['batch'] = 'seaad_dlpfc'
rna_sun.obs['batch'] = 'sun'

adata = ad.concat([rna_mtg, rna_dlpfc, rna_sun])
scvi.model.SCVI.setup_anndata(adata, layer="UMIs", batch_key="batch")

model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")

model.train()

SCVI_LATENT_KEY = "X_scVI"
adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()

adata.obs = adata.obs.astype(str)
adata.write_h5ad('/home/icb/zihe.zheng/projects/microglia/data/integrated_only_rna.h5ad')

model.save('/home/icb/zihe.zheng/projects/microglia/data/scvi_model/')






