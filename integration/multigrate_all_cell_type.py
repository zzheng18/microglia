import scanpy as sc
import multigrate as mtg
import anndata as ad

# # process atac_mtg
# n_top_genes = 5000

# atac_mtg = sc.read_h5ad('/lustre/groups/ml01/projects/2024_microglia_zihe.zheng/SEAAD/SEAAD_MTG_ATACseq_final-nuclei.2024-02-13.h5ad')
# # remove mislabeled obs
# atac_mtg = atac_mtg[~atac_mtg.obs['Supertype'].isin(['n_genes', 'fraction_mito', 'doublet_score','cluster_fraction_mito_flag','cluster_doublet_score_flag','cluster_fraction_ribo_flag', 'cluster_n_genes_flag'])]

# sc.pp.normalize_per_cell(atac_mtg, counts_per_cell_after=1e4)
# sc.pp.log1p(atac_mtg)
# sc.pp.highly_variable_genes(atac_mtg, n_top_genes=n_top_genes, subset=True)



rna_mtg = sc.read_h5ad('/lustre/groups/ml01/projects/2024_microglia_zihe.zheng/SEAAD/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13_count_hvg_aligned.h5ad')
atac_mtg = sc.read_h5ad('/lustre/groups/ml01/projects/2024_microglia_zihe.zheng/SEAAD/SEAAD_MTG_ATACseq_final-nuclei.2024-02-13_hvg_aligned.h5ad')
rna_dlpfc = sc.read_h5ad('/lustre/groups/ml01/projects/2024_microglia_zihe.zheng/SEAAD/SEAAD_DLPFC_RNAseq_final-nuclei.2024-02-13_count_hvg_aligned.h5ad')
rna_sun = sc.read_h5ad('/lustre/groups/ml01/projects/2024_microglia_zihe.zheng/rosmap/dlpfc2/raw_h5ad/all_hvg_seaad_aligned.h5ad')

print(rna_mtg)
print(rna_dlpfc)
print(rna_sun)
print(atac_mtg)


rna_mtg.obs['batch'] = 'seaad_mtg'
atac_mtg.obs['batch']= 'seaad_mtg'
rna_dlpfc.obs['batch'] = 'seaad_dlpfc'
rna_sun.obs['batch'] = 'rosmap'

rna_mtg.obs['Modality'] = 'rna'
atac_mtg.obs['Modality']= 'atac'
rna_dlpfc.obs['Modality'] = 'rna'
rna_sun.obs['Modality'] = 'rna'

rna_samp = rna_mtg.obs.sample_id.unique()
atac_samp = atac_mtg.obs.sample_id.unique()
overlap = set(rna_samp).intersection(atac_samp)

rna_mtg_overlap = rna_mtg[rna_mtg.obs.sample_id.isin(overlap)]
rna_mtg_overlap = rna_mtg_overlap[rna_mtg_overlap.obs_names.sort_values(), :]
rna_mtg_single = rna_mtg[~rna_mtg.obs.sample_id.isin(overlap)]

atac_mtg_overlap = atac_mtg[atac_mtg.obs.sample_id.isin(overlap)]
atac_mtg_overlap = atac_mtg_overlap[atac_mtg_overlap.obs_names.sort_values(), :]
atac_mtg_single = atac_mtg[~atac_mtg.obs.sample_id.isin(overlap)]

rna_mtg_overlap.obs['sub_batch'] = 'mtg_overlap'
rna_mtg_single.obs['sub_batch'] = 'mtg_single_rna'
atac_mtg_single.obs['sub_batch'] = 'mtg_single_atac'
rna_dlpfc.obs['sub_batch'] = 'dlpfc_rna'
rna_sun.obs['sub_batch'] = 'rosmap_rna'


adata = mtg.data.organize_multimodal_anndatas(
    adatas=[[rna_mtg_overlap, rna_mtg_single, rna_dlpfc, rna_sun], [atac_mtg_overlap, None, None, None]],
    layers=[['UMIs', 'UMIs','UMIs', None], [None, None, None, None]],
)

print(adata)

rna_indices_end = rna_mtg_overlap.shape[1]
mtg.model.MultiVAE.setup_anndata(
    adata,
    categorical_covariate_keys=['Modality', 'Donor_ID'],
    rna_indices_end=rna_indices_end,
)

vae = mtg.model.MultiVAE(
    adata,
    losses=["nb", "mse"],
    loss_coefs={
        "integ": 500,
    },
    integrate_on="Modality",
    alignment_type="marginal",
    modality_alignment="MMD",
)

vae.train(lr=1e-5, precision=32, max_epochs = 50)
vae.get_model_output()

vae.save('/home/icb/zihe.zheng/projects/microglia/model/seaad_rosmap_vae_model/')

adata.obs = adata.obs.astype(str)
sc.pp.neighbors(adata, use_rep="X_multigrate")
sc.tl.umap(adata)
adata.write_h5ad('/lustre/groups/ml01/projects/2024_microglia_zihe.zheng/integrated_data/seaad_rosmap_integrated.h5ad')