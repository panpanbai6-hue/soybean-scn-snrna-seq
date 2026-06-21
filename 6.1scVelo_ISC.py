###############################################################################
# Step 6.1: scVelo analysis for vascular and syncytium-related cells at 1 dpi
###############################################################################

import os
import scanpy as sc
import anndata as ad
from scipy import io
import numpy as np
import pandas as pd
import scvelo as scv

###############################################################################
# 1. Global settings
###############################################################################

np.random.seed(1)

scv.settings.verbosity = 3
scv.settings.set_figure_params(
    "scvelo",
    facecolor="white",
    dpi=300,
    frameon=False
)

transition_dir = "5_transition_objects"
loom_dir = "loom"
output_dir = "6_scvelo_vascular_syn_1dpi"

os.makedirs(output_dir, exist_ok=True)

cell_type_order = [
    "Pericycle",
    "Phloem",
    "Procambium",
    "Syncytium",
    "Xylem"
]

cell_type_palette = [
    "#556B2F",
    "#9ddcdc",
    "#C03B4C",
    "#456F97",
    "#ec610a"
]

###############################################################################
# 2. Forrest
###############################################################################

cultivar = "Forrest"

cultivar_dir = os.path.join(
    transition_dir,
    cultivar
)

cultivar_output_dir = os.path.join(
    output_dir,
    cultivar
)

os.makedirs(cultivar_output_dir, exist_ok=True)

scv.settings.figdir = cultivar_output_dir

###############################################################################
# 2.1 Load Seurat-exported matrix, metadata and embeddings
###############################################################################

X = io.mmread(
    os.path.join(
        cultivar_dir,
        "Forrest_vascular_syn_1dpi.mtx"
    )
)

adata = ad.AnnData(
    X=X.transpose().tocsr()
)

cell_meta = pd.read_csv(
    os.path.join(
        cultivar_dir,
        "Forrest_vascular_syn_1dpi_metadata.csv"
    )
)

cell_meta["barcode"] = cell_meta["barcode"].astype(str)
cell_meta.index = cell_meta["barcode"]

adata.obs = cell_meta

gene_names = pd.read_csv(
    os.path.join(
        cultivar_dir,
        "Forrest_vascular_syn_1dpi_gene_names.csv"
    ),
    header=None
).iloc[:, 0].astype(str).tolist()

gene_names = [
    gene.replace("soybean_", "", 1).replace("ann1.", "", 1)
    for gene in gene_names
]

adata.var_names = gene_names
adata.var_names_make_unique()

pca_file = os.path.join(
    cultivar_dir,
    "Forrest_vascular_syn_1dpi_pca.csv"
)

if os.path.exists(pca_file):
    pca = pd.read_csv(pca_file)
    if pca.shape[0] == adata.n_obs:
        adata.obsm["X_pca"] = pca.to_numpy()

umap_file = os.path.join(
    cultivar_dir,
    "Forrest_vascular_syn_1dpi_umap.csv"
)

if os.path.exists(umap_file):
    umap = pd.read_csv(umap_file)
    if umap.shape[0] == adata.n_obs:
        adata.obsm["X_umap"] = umap.iloc[:, 0:2].to_numpy()

if "X_umap" not in adata.obsm:
    adata.obsm["X_umap"] = np.vstack(
        (
            adata.obs["UMAP_1"].to_numpy(),
            adata.obs["UMAP_2"].to_numpy()
        )
    ).T

adata.obs["cell_type"] = pd.Categorical(
    adata.obs["cell_type"],
    categories=cell_type_order,
    ordered=True
)

###############################################################################
# 2.2 Load Forrest loom files
###############################################################################

F1_1 = sc.read(
    os.path.join(
        loom_dir,
        "F1_1.loom"
    ),
    cache=True
)

F1_2 = sc.read(
    os.path.join(
        loom_dir,
        "F1_2.loom"
    ),
    cache=True
)

F1_1 = F1_1[
    :,
    [
        not str(gene).startswith("scn") and
        not str(gene).startswith("SCN")
        for gene in F1_1.var_names
    ]
].copy()

F1_2 = F1_2[
    :,
    [
        not str(gene).startswith("scn") and
        not str(gene).startswith("SCN")
        for gene in F1_2.var_names
    ]
].copy()

F1_1.var_names = [
    str(gene).replace("soybean_", "", 1).replace("ann1.", "", 1)
    for gene in F1_1.var_names
]

F1_2.var_names = [
    str(gene).replace("soybean_", "", 1).replace("ann1.", "", 1)
    for gene in F1_2.var_names
]

F1_1.var_names_make_unique()
F1_2.var_names_make_unique()

F1_1.obs_names = [
    str(bc).split(":")[-1][:-1] + "_1"
    for bc in F1_1.obs_names
]

F1_2.obs_names = [
    str(bc).split(":")[-1][:-1] + "_2"
    for bc in F1_2.obs_names
]

F1_1.obs["library"] = "F1_1"
F1_2.obs["library"] = "F1_2"

ldata = ad.concat(
    [F1_1, F1_2],
    axis=0,
    join="outer",
    index_unique=None
)

ldata.var_names_make_unique()

print("Forrest Seurat cells:", adata.n_obs)
print("Forrest loom cells:", ldata.n_obs)
print("Forrest common cells:", len(adata.obs_names.intersection(ldata.obs_names)))
print("Forrest common genes:", len(adata.var_names.intersection(ldata.var_names)))

###############################################################################
# 2.3 Merge Seurat metadata with spliced and unspliced layers
###############################################################################

adata_all = scv.utils.merge(
    adata,
    ldata
)

adata_all.write(
    os.path.join(
        cultivar_output_dir,
        "Forrest_vascular_syn_1dpi_merged_raw.h5ad"
    )
)

sc.pl.umap(
    adata_all,
    color=["cell_type"],
    frameon=False,
    save=False
)

###############################################################################
# 2.4 Forrest scVelo analysis for all vascular and syncytium-related cells
###############################################################################

scv.pl.proportions(
    adata_all,
    groupby="cell_type",
    save="_Forrest_spliced_unspliced_proportions.pdf"
)

scv.pp.filter_and_normalize(
    adata_all,
    min_shared_counts=20,
    n_top_genes=3000
)

scv.pp.moments(
    adata_all,
    n_pcs=25,
    n_neighbors=25
)

scv.tl.recover_dynamics(
    adata_all,
    n_jobs=15
)

scv.tl.velocity(
    adata_all,
    mode="dynamical"
)

cells_to_flip = adata_all.obs["cell_type"].isin(
    ["Syncytium", "Procambium"]
)

adata_all.layers["velocity"][cells_to_flip.values, :] *= -1

scv.tl.velocity_graph(
    adata_all,
    n_jobs=15
)

scv.pl.velocity_embedding_stream(
    adata_all,
    basis="umap",
    color="cell_type",
    title="",
    palette=cell_type_palette,
    dpi=300,
    legend_loc="right margin",
    save="_Forrest_vascular_syn_embedding_stream.svg"
)

adata_all.write(
    os.path.join(
        cultivar_output_dir,
        "Forrest_vascular_syn_1dpi_scvelo_all.h5ad"
    )
)

###############################################################################
# 3. PI88788
###############################################################################

cultivar = "PI88788"

cultivar_dir = os.path.join(
    transition_dir,
    cultivar
)

cultivar_output_dir = os.path.join(
    output_dir,
    cultivar
)

os.makedirs(cultivar_output_dir, exist_ok=True)

scv.settings.figdir = cultivar_output_dir

###############################################################################
# 3.1 Load Seurat-exported matrix, metadata and embeddings
###############################################################################

X = io.mmread(
    os.path.join(
        cultivar_dir,
        "PI88788_vascular_syn_1dpi.mtx"
    )
)

adata = ad.AnnData(
    X=X.transpose().tocsr()
)

cell_meta = pd.read_csv(
    os.path.join(
        cultivar_dir,
        "PI88788_vascular_syn_1dpi_metadata.csv"
    )
)

cell_meta["barcode"] = cell_meta["barcode"].astype(str)
cell_meta.index = cell_meta["barcode"]

adata.obs = cell_meta

gene_names = pd.read_csv(
    os.path.join(
        cultivar_dir,
        "PI88788_vascular_syn_1dpi_gene_names.csv"
    ),
    header=None
).iloc[:, 0].astype(str).tolist()

gene_names = [
    gene.replace("soybean_", "", 1).replace("ann1.", "", 1)
    for gene in gene_names
]

adata.var_names = gene_names
adata.var_names_make_unique()

pca_file = os.path.join(
    cultivar_dir,
    "PI88788_vascular_syn_1dpi_pca.csv"
)

if os.path.exists(pca_file):
    pca = pd.read_csv(pca_file)
    if pca.shape[0] == adata.n_obs:
        adata.obsm["X_pca"] = pca.to_numpy()

umap_file = os.path.join(
    cultivar_dir,
    "PI88788_vascular_syn_1dpi_umap.csv"
)

if os.path.exists(umap_file):
    umap = pd.read_csv(umap_file)
    if umap.shape[0] == adata.n_obs:
        adata.obsm["X_umap"] = umap.iloc[:, 0:2].to_numpy()

if "X_umap" not in adata.obsm:
    adata.obsm["X_umap"] = np.vstack(
        (
            adata.obs["UMAP_1"].to_numpy(),
            adata.obs["UMAP_2"].to_numpy()
        )
    ).T

adata.obs["cell_type"] = pd.Categorical(
    adata.obs["cell_type"],
    categories=cell_type_order,
    ordered=True
)

###############################################################################
# 3.2 Load PI88788 loom files
###############################################################################

P1_1 = sc.read(
    os.path.join(
        loom_dir,
        "P1_1.loom"
    ),
    cache=True
)

P1_2 = sc.read(
    os.path.join(
        loom_dir,
        "P1_2.loom"
    ),
    cache=True
)

P1_1 = P1_1[
    :,
    [
        not str(gene).startswith("scn") and
        not str(gene).startswith("SCN")
        for gene in P1_1.var_names
    ]
].copy()

P1_2 = P1_2[
    :,
    [
        not str(gene).startswith("scn") and
        not str(gene).startswith("SCN")
        for gene in P1_2.var_names
    ]
].copy()

P1_1.var_names = [
    str(gene).replace("soybean_", "", 1).replace("ann1.", "", 1)
    for gene in P1_1.var_names
]

P1_2.var_names = [
    str(gene).replace("soybean_", "", 1).replace("ann1.", "", 1)
    for gene in P1_2.var_names
]

P1_1.var_names_make_unique()
P1_2.var_names_make_unique()

P1_1.obs_names = [
    str(bc).split(":")[-1][:-1] + "_1"
    for bc in P1_1.obs_names
]

P1_2.obs_names = [
    str(bc).split(":")[-1][:-1] + "_2"
    for bc in P1_2.obs_names
]

P1_1.obs["library"] = "P1_1"
P1_2.obs["library"] = "P1_2"

ldata = ad.concat(
    [P1_1, P1_2],
    axis=0,
    join="outer",
    index_unique=None
)

ldata.var_names_make_unique()

print("PI88788 Seurat cells:", adata.n_obs)
print("PI88788 loom cells:", ldata.n_obs)
print("PI88788 common cells:", len(adata.obs_names.intersection(ldata.obs_names)))
print("PI88788 common genes:", len(adata.var_names.intersection(ldata.var_names)))

###############################################################################
# 3.3 Merge Seurat metadata with spliced and unspliced layers
###############################################################################

adata_all = scv.utils.merge(
    adata,
    ldata
)

adata_all.write(
    os.path.join(
        cultivar_output_dir,
        "PI88788_vascular_syn_1dpi_merged_raw.h5ad"
    )
)

sc.pl.umap(
    adata_all,
    color=["cell_type"],
    frameon=False,
    save=False
)

###############################################################################
# 3.4 PI88788 scVelo analysis for all vascular and syncytium-related cells
###############################################################################

scv.pl.proportions(
    adata_all,
    groupby="cell_type",
    save="_PI88788_spliced_unspliced_proportions.pdf"
)

scv.pp.filter_and_normalize(
    adata_all,
    min_shared_counts=20,
    n_top_genes=3000
)

scv.pp.moments(
    adata_all,
    n_pcs=25,
    n_neighbors=25
)

scv.tl.recover_dynamics(
    adata_all,
    n_jobs=15
)

scv.tl.velocity(
    adata_all,
    mode="dynamical"
)

cells_to_flip = adata_all.obs["cell_type"].isin(
    ["Syncytium", "Procambium"]
)

adata_all.layers["velocity"][cells_to_flip.values, :] *= -1

scv.tl.velocity_graph(
    adata_all,
    n_jobs=15
)

scv.pl.velocity_embedding_stream(
    adata_all,
    basis="umap",
    color="cell_type",
    title="",
    palette=cell_type_palette,
    dpi=300,
    legend_loc="right margin",
    save="_PI88788_vascular_syn_embedding_stream.svg"
)

adata_all.write(
    os.path.join(
        cultivar_output_dir,
        "PI88788_vascular_syn_1dpi_scvelo_all.h5ad"
    )
)


###############################################################################
# 4. Williams82
###############################################################################

cultivar = "Williams82"

cultivar_dir = os.path.join(
    transition_dir,
    cultivar
)

cultivar_output_dir = os.path.join(
    output_dir,
    cultivar
)

os.makedirs(cultivar_output_dir, exist_ok=True)

scv.settings.figdir = cultivar_output_dir

###############################################################################
# 4.1 Load Seurat-exported matrix, metadata and embeddings
###############################################################################

X = io.mmread(
    os.path.join(
        cultivar_dir,
        "Williams82_vascular_syn_1dpi.mtx"
    )
)

adata = ad.AnnData(
    X=X.transpose().tocsr()
)

cell_meta = pd.read_csv(
    os.path.join(
        cultivar_dir,
        "Williams82_vascular_syn_1dpi_metadata.csv"
    )
)

cell_meta["barcode"] = cell_meta["barcode"].astype(str)
cell_meta.index = cell_meta["barcode"]

adata.obs = cell_meta

gene_names = pd.read_csv(
    os.path.join(
        cultivar_dir,
        "Williams82_vascular_syn_1dpi_gene_names.csv"
    ),
    header=None
).iloc[:, 0].astype(str).tolist()

gene_names = [
    gene.replace("soybean_", "", 1).replace("ann1.", "", 1)
    for gene in gene_names
]

adata.var_names = gene_names
adata.var_names_make_unique()

pca_file = os.path.join(
    cultivar_dir,
    "Williams82_vascular_syn_1dpi_pca.csv"
)

if os.path.exists(pca_file):
    pca = pd.read_csv(pca_file)
    if pca.shape[0] == adata.n_obs:
        adata.obsm["X_pca"] = pca.to_numpy()

umap_file = os.path.join(
    cultivar_dir,
    "Williams82_vascular_syn_1dpi_umap.csv"
)

if os.path.exists(umap_file):
    umap = pd.read_csv(umap_file)
    if umap.shape[0] == adata.n_obs:
        adata.obsm["X_umap"] = umap.iloc[:, 0:2].to_numpy()

if "X_umap" not in adata.obsm:
    adata.obsm["X_umap"] = np.vstack(
        (
            adata.obs["UMAP_1"].to_numpy(),
            adata.obs["UMAP_2"].to_numpy()
        )
    ).T

adata.obs["cell_type"] = pd.Categorical(
    adata.obs["cell_type"],
    categories=cell_type_order,
    ordered=True
)

###############################################################################
# 4.2 Load Williams82 loom files
###############################################################################

W1_1 = sc.read(
    os.path.join(
        loom_dir,
        "W1_1.loom"
    ),
    cache=True
)

W1_2 = sc.read(
    os.path.join(
        loom_dir,
        "W1_2.loom"
    ),
    cache=True
)

W1_1 = W1_1[
    :,
    [
        not str(gene).startswith("scn") and
        not str(gene).startswith("SCN")
        for gene in W1_1.var_names
    ]
].copy()

W1_2 = W1_2[
    :,
    [
        not str(gene).startswith("scn") and
        not str(gene).startswith("SCN")
        for gene in W1_2.var_names
    ]
].copy()

W1_1.var_names = [
    str(gene).replace("soybean_", "", 1).replace("ann1.", "", 1)
    for gene in W1_1.var_names
]

W1_2.var_names = [
    str(gene).replace("soybean_", "", 1).replace("ann1.", "", 1)
    for gene in W1_2.var_names
]

W1_1.var_names_make_unique()
W1_2.var_names_make_unique()

W1_1.obs_names = [
    str(bc).split(":")[-1][:-1] + "_1"
    for bc in W1_1.obs_names
]

W1_2.obs_names = [
    str(bc).split(":")[-1][:-1] + "_2"
    for bc in W1_2.obs_names
]

W1_1.obs["library"] = "W1_1"
W1_2.obs["library"] = "W1_2"

ldata = ad.concat(
    [W1_1, W1_2],
    axis=0,
    join="outer",
    index_unique=None
)

ldata.var_names_make_unique()

print("Williams82 Seurat cells:", adata.n_obs)
print("Williams82 loom cells:", ldata.n_obs)
print("Williams82 common cells:", len(adata.obs_names.intersection(ldata.obs_names)))
print("Williams82 common genes:", len(adata.var_names.intersection(ldata.var_names)))

###############################################################################
# 4.3 Merge Seurat metadata with spliced and unspliced layers
###############################################################################

adata_all = scv.utils.merge(
    adata,
    ldata
)

adata_all.write(
    os.path.join(
        cultivar_output_dir,
        "Williams82_vascular_syn_1dpi_merged_raw.h5ad"
    )
)

sc.pl.umap(
    adata_all,
    color=["cell_type"],
    frameon=False,
    save=False
)

###############################################################################
# 4.4 Williams82 scVelo analysis for all vascular and syncytium-related cells
###############################################################################

scv.pl.proportions(
    adata_all,
    groupby="cell_type",
    save="_Williams82_spliced_unspliced_proportions.pdf"
)

scv.pp.filter_and_normalize(
    adata_all,
    min_shared_counts=20,
    n_top_genes=3000
)

scv.pp.moments(
    adata_all,
    n_pcs=30,
    n_neighbors=20
)

scv.tl.recover_dynamics(
    adata_all,
    n_jobs=15
)

scv.tl.velocity(
    adata_all,
    mode="dynamical"
)

cells_to_flip = adata_all.obs["cell_type"].isin(
    ["Syncytium", "Procambium"]
)

adata_all.layers["velocity"][cells_to_flip.values, :] *= -1

scv.tl.velocity_graph(
    adata_all,
    n_jobs=15
)

scv.pl.velocity_embedding_stream(
    adata_all,
    basis="umap",
    color="cell_type",
    title="",
    palette=cell_type_palette,
    dpi=300,
    legend_loc="right margin",
    save="_Williams82_vascular_syn_embedding_stream.svg"
)

adata_all.write(
    os.path.join(
        cultivar_output_dir,
        "Williams82_vascular_syn_1dpi_scvelo_all.h5ad"
    )
)