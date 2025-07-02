import matplotlib.pyplot as plt
import scanpy as sc
import scvi
from scvi.external import SysVI
import anndata as ad
import pandas as pd

scvi.settings.seed = 0

adata = sc.read_h5ad('/lustre/groups/ml01/projects/2024_microglia_zihe.zheng/integrated_data/human_mouse_not_integrated.h5ad')

## integration
# Setup adata for training
SysVI.setup_anndata(
    adata=adata,
    batch_key="system",
    categorical_covariate_keys=["batch"],
)

# Example showing how to turn on embedding of all categorical covariates
model = SysVI(
    adata=adata,
    embed_categorical_covariates=True,
)

# Train
max_epochs = 200
model.train(
    max_epochs=max_epochs, check_val_every_n_epoch=1, plan_kwargs={"kl_weight": 5, "z_distance_cycle_weight": 2}
)

# Make detailed plot after N epochs
epochs_detail_plot = 100

# Losses to plot
losses = [
    "reconstruction_loss_train",
    "kl_local_train",
    "cycle_loss_train",
]
fig, axs = plt.subplots(2, len(losses), figsize=(len(losses) * 3, 4))
for ax_i, l_train in enumerate(losses):
    l_val = l_train.replace("_train", "_validation")
    l_name = l_train.replace("_train", "")
    # Change idx of epochs to start with 1
    l_val_values = model.trainer.logger.history[l_val].copy()
    l_val_values.index = l_val_values.index + 1
    l_train_values = model.trainer.logger.history[l_train].copy()
    l_train_values.index = l_train_values.index + 1
    for l_values, c, alpha, dp in [
        (l_train_values, "tab:blue", 1, epochs_detail_plot),
        (l_val_values, "tab:orange", 0.5, epochs_detail_plot),
    ]:
        axs[0, ax_i].plot(l_values.index, l_values.values.ravel(), c=c, alpha=alpha)
        axs[0, ax_i].set_title(l_name)
        axs[1, ax_i].plot(l_values.index[dp:], l_values.values.ravel()[dp:], c=c, alpha=alpha)

fig.tight_layout()
fig.savefig("/home/icb/zihe.zheng/projects/microglia/plots/integration_sysvi_all.pdf")

# Get embedding - save it into X of new AnnData
embed = model.get_latent_representation(adata=adata)
embed = sc.AnnData(embed, obs=adata.obs)
# Make system categorical for plotting below
embed.obs["system"] = embed.obs["system"]
# Compute UMAP
sc.pp.neighbors(embed, use_rep="X")
sc.tl.umap(embed)

# embed.obs['donor'] = embed.obs['donor'].astype(str)
embed.write_h5ad('/lustre/groups/ml01/projects/2024_microglia_zihe.zheng/integrated_data/human_mouse.h5ad')

# Obs columns to color by
cols = ["system", "state", 'fine_cluster', 'paper_cluster']

# One plot per obs column used for coloring
fig, axs = plt.subplots(len(cols), 1, figsize=(3, 3 * len(cols)))
for col, ax in zip(cols, axs, strict=False):
    sc.pl.embedding(
        embed,
        "X_umap",
        color=col,
        s=10,
        ax=ax,
        show=False,
        sort_order=False,
        frameon=False,
    )
fig.tight_layout()
fig.savefig("/home/icb/zihe.zheng/projects/microglia/plots/integration_sysvi_umap_all.pdf", dpi=300)