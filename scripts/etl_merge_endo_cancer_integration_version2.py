use_gpu = True # whether to use GPU acceleration for the SCVI data integration
# Remember to give a folder name to the run in line 48

"""
Merges the data of endometriosis and cancer datasets, saves merged data as AnnData.
A dataframe with the target (Cancer/EMS/Normal) along with gene counts/metadata is saved as parquet.
A dataframe with the target (Cancer/EMS/Normal) is saved as parquet.

Source data tar files need to be saved in the following locations:
    - data/raw/cancer/GSE184880_RAW.tar 
    - data/raw/endo/GSE214411_RAW.tar
    
The processed files get stored in a folder with this format: data/processed/name
Rename file in line 48 before running code
"""

# Import libraries
import os
import tarfile
import shutil 
import time

import scanpy as sc
import anndata as ad
import scvi

from sklearn.model_selection import StratifiedGroupKFold

# Define function for extracting highly variable features 
def reduce_features(adata):
    # Assign a copy
    adata_copy = adata.copy()
    # Normalize
    sc.pp.normalize_total(adata_copy, target_sum=1e4)
    # Logarithmize
    sc.pp.log1p(adata_copy)
    # Identify highly variable genes
    sc.pp.highly_variable_genes(adata_copy, min_mean=0.0125, max_mean=3, min_disp=0.5)
    # Plug the highly variable genes back into the original adata
    adata = adata[:, adata_copy.var.highly_variable]
    return adata

print(os.path.join(os.path.dirname(__file__), ".."))

project_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), ".."))

# Create folders for processed data
processed_data_path = os.path.normpath(os.path.join(project_dir, "data", "processed", "v3_integrated"))
os.makedirs(processed_data_path, exist_ok=False) # exist_ok = false prevents overwriting existing data folder


# Source folders
dir_endo    = os.path.join(project_dir, "data", "raw", "endo") # path for directory containing endometriosis data .tar file.
dir_cancer  = os.path.join(project_dir, "data", "raw", "cancer") # path for directory containing cancer data .tar file.


## Endometriosis Data
# Extract files to a directory
data_dir = os.path.join(dir_endo, "temp")
with tarfile.open(os.path.join(dir_endo, "GSE214411_RAW.tar"), "r") as tar:
    tar.extractall(path=data_dir)
# List all prefixes (one per sample)
prefixes = [
    "GSM6605431_EMS1_",
    "GSM6605432_EMS2_",
    "GSM6605433_EMS3_",
    "GSM6605434_EMS4_",
    "GSM6605435_EMS5_",
    "GSM6605436_EMS6_",
    "GSM6605437_N1_",
    "GSM6605438_N2_",
    "GSM6605439_N3_",
    "GSM6605440_N4_",
    "GSM7277296_N-5_",
    "GSM7277297_N-6_",
    "GSM7277298_N-7_",
]

adatas = []
for prefix in prefixes:
    adata = sc.read_10x_mtx(
        data_dir,
        prefix=prefix,
        var_names="gene_symbols",
        make_unique=True,
        cache=True,
    )
    adata.obs["sample"] = prefix.rstrip("_")
    adatas.append(adata)
    
# Delete extracted files
shutil.rmtree(data_dir)

# Concatenate all Objects inside adatas into one
adata_endo = ad.concat(
    adatas, label="batch", index_unique="-", keys=[p.rstrip("_") for p in prefixes]
)

# Basic filtering: number of genes/cell and number of cells/samples
sc.pp.filter_cells(adata_endo, min_genes=200)
sc.pp.filter_genes(adata_endo, min_cells=3)
# annotate the group of mitochondrial genes as "mt"
adata_endo.var["mt"] = adata_endo.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata_endo, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)
# The actual filtering (refer to ETL notebook for details)
adata_endo = adata_endo[
    adata_endo.obs.n_genes_by_counts < 2500, :
]  # to remove doublets (two or more cells in one droplet)
adata_endo = adata_endo[
    adata_endo.obs.pct_counts_mt < 20, :
].copy()  # to remove cell with high mitochondrial activity (dying  cells)

## Ovarian Cancer Data
# For loop to create an adata for the cancer dataset
# Extract files to a directory

data_dir = os.path.join(dir_cancer, "temp")
with tarfile.open(os.path.join(dir_cancer, "GSE184880_RAW.tar"), "r") as tar:
    tar.extractall(path=data_dir)
# List all prefixes (one per sample)
prefixes = [
    "GSM5599220_Norm1",
    "GSM5599221_Norm2",
    "GSM5599222_Norm3",
    "GSM5599223_Norm4",
    "GSM5599224_Norm5",
    "GSM5599225_Cancer1",
    "GSM5599226_Cancer2",
    "GSM5599227_Cancer3",
    "GSM5599228_Cancer4",
    "GSM5599229_Cancer5",
    "GSM5599230_Cancer6",
    "GSM5599231_Cancer7",
]

adatas = []
for prefix in prefixes:
    # Rename files to match Scanpy's expectations

    os.rename(
        os.path.join(data_dir, f"{prefix}.matrix.mtx.gz"),
        os.path.join(data_dir, f"{prefix}_matrix.mtx.gz"),
    )
    os.rename(
        os.path.join(data_dir, f"{prefix}.barcodes.tsv.gz"),
        os.path.join(data_dir, f"{prefix}_barcodes.tsv.gz"),
    )
    os.rename(
        os.path.join(data_dir, f"{prefix}.genes.tsv.gz"),
        os.path.join(data_dir, f"{prefix}_features.tsv.gz"),
    )
    prefix = prefix + "_"  # add underscore
    adata = sc.read_10x_mtx(
        data_dir,
        prefix=prefix,
        var_names="gene_symbols",
        make_unique=True,
        cache=True,
    )
    adata.obs["sample"] = prefix.rstrip("_")
    adatas.append(adata)
        
# Delete extracted files
shutil.rmtree(data_dir)

# Concatenate all Objects inside adatas into one
adata_cancer = ad.concat(
    adatas, label="batch", index_unique="-", keys=[p.rstrip("_") for p in prefixes]
)

# Basic filtering: number of genes/cell and number of cells/samples
sc.pp.filter_cells(adata_cancer, min_genes=200)
sc.pp.filter_genes(adata_cancer, min_cells=3)
# annotate the group of mitochondrial genes as "mt"
adata_cancer.var["mt"] = adata_cancer.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata_cancer, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)
# The actual filtering, refer to ETL notebook for details
adata_cancer = adata_cancer[
    adata_cancer.obs.n_genes_by_counts < 3000, :
]  # to remove doublets (two or more cells in one droplet)
adata_cancer = adata_cancer[
    adata_cancer.obs.pct_counts_mt < 20, :
].copy()  # to remove cells with high mitochondrial activity (dying cells)

# Highly variable genes
adata_cancer = reduce_features(adata_cancer)
adata_endo = reduce_features(adata_endo)

## Merge the two datasets
adata_merged = ad.concat(
    [adata_endo, adata_cancer],
    axis=0,  # Concatenate along cells
    join="inner",  # Keep intersected genes only
    keys=["Endometriosis", "Cancer"],
    merge="unique",  # Handle overlapping metadata uniquely
    index_unique="-",  # Avoid duplicate observation names
    fill_value=0,  # Fill missing values wth zeros
)

# Assign datasets as a column
adata_merged.obs['dataset'] = adata_merged.obs.index.str.split('-').str[-1]
# Assign samples as a column
adata_merged.obs['samples'] = adata_merged.obs["sample"].str.split('-').str[0]

# Create target variable in obs
adata_merged.obs["target"] = (
    adata_merged.obs["sample"]
    .str.split("_").str[-1]
    .str[0]
    .map({"C": "Cancer", "E": "EMS", "N": "Normal"})
)

# Extract group and target information from adata.obs
groups = adata_merged.obs['samples']
y = adata_merged.obs['target']

# Create the splitter (n_splits=2 for train/test)
sgkf = StratifiedGroupKFold(n_splits=4, shuffle=True, random_state=42)

# Get the first split (train/test indices)
train_idx, test_idx = next(sgkf.split(adata_merged.X, y, groups))

# Subset the AnnData object for train and test
adata_train = adata_merged[train_idx].copy()
adata_test = adata_merged[test_idx].copy()

adata_train.layers["counts"] = adata_train.X
adata_test.layers["counts"] = adata_test.X

## Integrate to regress technical variations

# Setup model, fit, and normalise train dataset 
scvi.model.SCVI.setup_anndata(adata_train, layer="counts", 
                              categorical_covariate_keys= ["sample", "dataset"], # choosing datasets as a covariate
                              continuous_covariate_keys= ["pct_counts_mt", "total_counts"])
# Fit on train set only
model = scvi.model.SCVI(adata_train, n_layers=2, n_latent=30, gene_likelihood="nb")
if use_gpu == True:
    model.train(accelerator="gpu") 
else: model.train() 

# Get normalized counts for train set
norm_train = model.get_normalized_expression(
    adata_train, 
    library_size=1e4,
    return_numpy=False
)
# Insert integrated data to Anndata
latent_train = model.get_latent_representation()
adata_train.obsm["X_scVI"] = latent_train

# Setup model, fit, and normalise test dataset 
scvi.model.SCVI.setup_anndata(adata_test, layer="counts", 
                              categorical_covariate_keys= ["sample", "dataset"], # choosing datasets as a covariate
                              continuous_covariate_keys= ["pct_counts_mt", "total_counts"])
# Fit on train set only
model = scvi.model.SCVI(adata_test, n_layers=2, n_latent=30, gene_likelihood="nb")
if use_gpu == True:
    model.train(accelerator="gpu") 
else: model.train() 

# Get normalized counts for test set
norm_test = model.get_normalized_expression(
    adata_test,
    library_size=1e4, 
    return_numpy=False
)
# Insert integrated data to Anndata
latent_test = model.get_latent_representation()
adata_test.obsm["X_scVI"] = latent_test

# Create DataFrames with targets
norm_train_df = norm_train.join(adata_train.obs["target"])
norm_test_df = norm_test.join(adata_test.obs["target"])

# Save normalized datasets
norm_train_df.to_parquet(os.path.join(processed_data_path, "scvi_normalized_train.parquet"))
norm_test_df.to_parquet(os.path.join(processed_data_path, "scvi_normalized_test.parquet"))

# Merge data again 
adata_integrated = ad.concat(
    [adata_train, adata_test],
    axis=0,          # Concatenate along cells 
    join="outer",    # Keep all genes 
    merge="unique",  # Handle overlapping metadata uniquely
    fill_value=0     # Fill missing values wth zeros
)

# Save Anndata in case of visualisation ## cannot generally save since the latent space has different dimensions.
ad.AnnData.write_h5ad(adata_integrated, filename=os.path.join(processed_data_path, "anndata_integrated.h5ad"), compression='gzip')
