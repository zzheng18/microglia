import numpy as np
import scanpy as sc
import pandas as pd
import anndata as ad
import scvi
import matplotlib as plt
import torch

rosmap = sc.read_h5ad('/lustre/groups/ml01/projects/2024_microglia_zihe.zheng/rosmap/dlpfc2/raw_h5ad/all_hvg_rosmap_align.h5ad')
diverse_cohort = sc.read_h5ad('/lustre/groups/ml01/projects/2024_microglia_zihe.zheng/diverse_cohort/diverse_cohort_hvg_rosmap_annot_align.h5ad')
multi_region = sc.read_h5ad('/lustre/groups/ml01/projects/2024_microglia_zihe.zheng/rosmap_multiome/all_brain_regions_filt_preprocessed_scanpy_fullmatrix_annot_hvg_rosmap_annot_align.h5ad')
aging = sc.read_h5ad('/lustre/groups/ml01/projects/2024_microglia_zihe.zheng/rosmap_aging/PFC427_raw_data_hvg_rosmap_annot_align.h5ad')

rosmap.X = rosmap.X.astype(np.float32)
diverse_cohort.X = diverse_cohort.X.astype(np.float32)

adata = ad.concat([rosmap, diverse_cohort, multi_region, aging], join = 'outer')
adata.obs_names_make_unique()

scvi.model.SCVI.setup_anndata(adata, batch_key="dataset_batch")

model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")

model.train(max_epochs=50)

SCVI_LATENT_KEY = "X_scVI"
adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()

adata.obs = adata.obs.astype(str)

sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY)
sc.tl.umap(adata)

sc.pl.umap(adata, color = ['dataset_batch', 'aligned_cell_type'], show = False, save='_rosmap_all.pdf')

adata.write_h5ad('/lustre/groups/ml01/projects/2024_microglia_zihe.zheng/integrated_data/rosmap_everything_integrated.h5ad')

# model.save('/home/icb/zihe.zheng/projects/microglia/data/scvi_model/')

history = model.history

# Plot
plt.figure(figsize=(8, 5))

# Plot available metrics if they exist
if 'train_loss_epoch' in history:
    plt.plot(history['train_loss_epoch'], label='Training Loss')
if 'elbo_train' in history:
    plt.plot(history['elbo_train'], label='Train ELBO')
if 'reconstruction_loss_train' in history:
    plt.plot(history['reconstruction_loss_train'], label='train reconstruction loss')
if 'elbo_validation' in history:
    plt.plot(history['elbo_validation'], label='Validation ELBO')
if 'reconstruction_loss_validation' in history:
    plt.plot(history['reconstruction_loss_validation'], label='validation reconstruction loss')


plt.xlabel('Epoch')
plt.ylabel('Loss / ELBO')
plt.title('SCVI Training Metrics')
plt.legend()
plt.tight_layout()

# Save figure
plt.savefig('/home/icb/zihe.zheng/projects/microglia/plots/scvi_rosmap_losses.pdf', dpi=300)

