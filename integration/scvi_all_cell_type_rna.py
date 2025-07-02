import numpy as np
import scanpy as sc
import anndata as ad
import scvi
import matplotlib as plt
import torch

rna_mtg = sc.read_h5ad('/lustre/groups/ml01/projects/2024_microglia_zihe.zheng/SEAAD/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13_count_hvg_aligned.h5ad')
rna_dlpfc = sc.read_h5ad('/lustre/groups/ml01/projects/2024_microglia_zihe.zheng/SEAAD/SEAAD_DLPFC_RNAseq_final-nuclei.2024-02-13_count_hvg_aligned.h5ad')
rna_sun = sc.read_h5ad('/lustre/groups/ml01/projects/2024_microglia_zihe.zheng/rosmap/dlpfc2/raw_h5ad/all_hvg_seaad_aligned.h5ad')

# rna_sun.layers['UMIs'] = rna_sun.X.copy()

rna_mtg.obs['batch'] = 'seaad_mtg'
rna_dlpfc.obs['batch'] = 'seaad_dlpfc'
rna_sun.obs['batch'] = 'sun'

adata = ad.concat([rna_mtg, rna_dlpfc, rna_sun])
scvi.model.SCVI.setup_anndata(adata, batch_key="batch")

sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

sc.pl.umap(adata, color = ['batch'], show = False, save = '_seaad_rosmap_integrated_rna_before_int.pdf')
# plt.savefig('/home/icb/zihe.zheng/projects/microglia/plots/seaad_rosmap_integrated_rna_before_int.pdf')
# plt.close()

model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")

model.train()

model.plot_losses(save = '/home/icb/zihe.zheng/projects/microglia/plots/seaad_rosmap_integrated_rna_before_int.pdf')

SCVI_LATENT_KEY = "X_scVI"
adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()

adata.obs = adata.obs.astype(str)

sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata)

# align cell type labels
subclass_to_broad = {
    'L5 IT': 'Excitatory',
    'L6 IT': 'Excitatory',
    'Pvalb': 'Inhibitory',
    'Microglia-PVM': 'Microglia',
    'Astrocyte': 'Astrocytes',
    'L4 IT': 'Excitatory',
    'L2/3 IT': 'Excitatory',
    'L5/6 NP': 'Excitatory',
    'Sst': 'Inhibitory',
    'Vip': 'Inhibitory',
    'L6 IT Car3': 'Excitatory',
    'Lamp5': 'Inhibitory',
    'Oligodendrocyte': 'Oligodendrocytes',
    'Lamp5 Lhx6': 'Inhibitory',
    'L6 CT': 'Excitatory',
    'Chandelier': 'Inhibitory',
    'L6b': 'Excitatory',
    'OPC': 'OPCs',
    'Sncg': 'Inhibitory',
    'Endothelial': 'Endothelial',
    'Pax6': 'Inhibitory',  
    'VLMC': 'Endothelial',  
    'Sst Chodl': 'Inhibitory',
    'L5 ET': 'Excitatory'
}

# Map to new column
adata.obs['Subclass_aligned'] = adata.obs['Subclass'].map(subclass_to_broad)
adata.obs["subset"] = adata.obs["subset"].cat.rename_categories({"CUX2+": "Excitatory"})
adata.obs['Subclass_aligned'] = adata.obs['Subclass_aligned'].fillna(adata.obs['subset'])

sc.pl.umap(adata, color = ['batch', 'Subclass_aligned'], show = False, save='_seaad_rosmap_integrated_rna_after_int.pdf')

adata.write_h5ad('/home/icb/zihe.zheng/projects/microglia/data/seaad_rosmap_integrated_rna_scvi.h5ad')

model.save('/home/icb/zihe.zheng/projects/microglia/data/scvi_model/')



