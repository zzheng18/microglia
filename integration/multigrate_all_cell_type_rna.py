import scanpy as sc
import multigrate as mtg
import anndata as ad


rna_mtg = sc.read_h5ad('/lustre/groups/ml01/projects/2024_microglia_zihe.zheng/SEAAD/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13_count_hvg_aligned.h5ad')
rna_dlpfc = sc.read_h5ad('/lustre/groups/ml01/projects/2024_microglia_zihe.zheng/SEAAD/SEAAD_DLPFC_RNAseq_final-nuclei.2024-02-13_count_hvg_aligned.h5ad')
rna_sun = sc.read_h5ad('/lustre/groups/ml01/projects/2024_microglia_zihe.zheng/rosmap/dlpfc2/raw_h5ad/all_hvg_seaad_aligned.h5ad')

print(rna_mtg)
print(rna_dlpfc)
print(rna_sun)

rna_mtg.obs['batch'] = 'seaad_mtg'
rna_dlpfc.obs['batch'] = 'seaad_dlpfc'
rna_sun.obs['batch'] = 'rosmap'

rna_mtg.obs['Modality'] = 'rna'
rna_dlpfc.obs['Modality'] = 'rna'
rna_sun.obs['Modality'] = 'rna'


adata = mtg.data.organize_multimodal_anndatas(
    adatas=[[rna_mtg, rna_dlpfc, rna_sun]],
    layers=[['UMIs','UMIs', None]],
)

print(adata)

rna_indices_end = rna_mtg.shape[1] + rna_dlpfc.shape[1] + rna_sun.shape[1]
mtg.model.MultiVAE.setup_anndata(
    adata,
    categorical_covariate_keys=['batch'],
    rna_indices_end=rna_indices_end,
)

vae = mtg.model.MultiVAE(
    adata,
    losses=["nb"]
)

vae.train(lr=1e-5, max_epochs = 50)
vae.get_model_output()

vae.save('/home/icb/zihe.zheng/projects/microglia/model/seaad_rosmap_vae_model_rna/')

vae.plot_losses(save = '/home/icb/zihe.zheng/projects/microglia/plots/integration_rna.pdf')

adata.obs = adata.obs.astype(str)
sc.pp.neighbors(adata, use_rep="X_multigrate")
sc.tl.umap(adata)
adata.write_h5ad('/lustre/groups/ml01/projects/2024_microglia_zihe.zheng/integrated_data/seaad_rosmap_integrated_rna.h5ad')