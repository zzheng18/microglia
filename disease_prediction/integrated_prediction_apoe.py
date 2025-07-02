import scanpy as sc
import multimil as mtm
import anndata as ad
import numpy as np
import random
import pandas as pd
import matplotlib.pyplot as plt
import scvi
from sklearn.metrics import classification_report, accuracy_score, f1_score, confusion_matrix, recall_score, balanced_accuracy_score
import scipy.stats as st

# read integrated data
adata = sc.read_h5ad('/lustre/groups/ml01/projects/2024_microglia_zihe.zheng/integrated_data/seaad_rosmap_integrated_rna_harmonized.h5ad')

# define parameters
sample_key = 'Donor_ID'
classification_keys = "apoe_label"
z_dim = 16
categorical_covariate_keys = [classification_keys] + [sample_key]
seed = 0

# create new adata with multigrate embeddings
new_adata = ad.AnnData(X=adata.obsm['X_multigrate'])
new_adata.obs = adata.obs.copy()
new_adata.obsm = adata.obsm.copy()
new_adata = new_adata[~new_adata.obs.apoe_label.isna()]

# stratified sampling by category 
query_samples = []
query_proportion = 0.2
for cat in new_adata.obs[classification_keys].unique():
    adata_cat = new_adata[new_adata.obs[classification_keys] == cat]
    samples = list(np.unique(adata_cat.obs[sample_key]))
    n_samples = len(samples)
    random.seed(seed)
    query_samples.extend(random.sample(samples, int(n_samples * query_proportion)))
    
print(query_samples)

# separate query and training adata
query = new_adata[new_adata.obs[sample_key].isin(query_samples)].copy()
adata = new_adata[~new_adata.obs[sample_key].isin(query_samples)].copy()
query.obs["ref"] = "query"
adata.obs["ref"] = "reference"

# sort index in both training adata and query 
idx = adata.obs[sample_key].sort_values().index
adata = adata[idx].copy()
idx = query.obs[sample_key].sort_values().index
query = query[idx].copy()

# train mil model on training set
scvi.settings.seed = seed
mtm.model.MILClassifier.setup_anndata(
        adata,
        categorical_covariate_keys=categorical_covariate_keys,
    )
print('seeeeeeeeeeeeeeeeeeeed', seed)
mil = mtm.model.MILClassifier(
    adata,
    classification=[classification_keys],
    z_dim=z_dim,
    sample_key=sample_key,
    # class_loss_coef = 0.1,
    # dropout = 0.05,
    # n_layers_cell_aggregator = 16,
    # n_layers_classifier = 4,
)

mil.train(lr=1e-3, max_epochs = 500, progress_bar_refresh_rate=0.5, train_size = 0.8, batch_size = 256, check_val_every_n_epoch=1)
mil.plot_losses(save = '/home/icb/zihe.zheng/projects/microglia/plots/seaad_rosmap_apoe_mil_rna.pdf')
mil.get_model_output()

mil.save('/home/icb/zihe.zheng/projects/microglia/model_mil/seaad_rosmap_apoe_mil_rna/')

# set up query data and infernece on query data
mtm.model.MILClassifier.setup_anndata(
    query,
    categorical_covariate_keys=categorical_covariate_keys,
)
new_model = mtm.model.MILClassifier.load_query_data(query, mil)
new_model.get_model_output()

# attach query to adata and save
adata_both = ad.concat([adata, query]) 
adata_both.write_h5ad('/lustre/groups/ml01/projects/2024_microglia_zihe.zheng/integrated_data/seaad_rosmap_integrated_rna_harmonized_apoe.h5ad')

# calculate scores on the query dataset
accuracy = accuracy_score(query.obs["apoe_label"], query.obs["predicted_apoe_label"])
print('accuracy', accuracy)

score = f1_score(query.obs["apoe_label"], query.obs["predicted_apoe_label"], average='macro')
print('F1_score', score)

cm = confusion_matrix(query.obs["apoe_label"], query.obs["predicted_apoe_label"])
print(cm)

balanced_acc = balanced_accuracy_score(query.obs["apoe_label"], query.obs["predicted_apoe_label"])
print('balanced accuracy', balanced_acc)

recalls = recall_score(query.obs["apoe_label"], query.obs["predicted_apoe_label"], average=None)
g_mean = np.prod(recalls) ** (1.0 / len(recalls))
print('g_mean', g_mean)

