"""Returns the merged data of endometriosis and cancer datasets.

Args:
    dir_endo (_type_): path for directory containing endometriosis data .tar files.
    dir_cancer (_type_): path for directory containing cancer data .tar files.

Returns:
    Anndata object: merged data of both datasets.
"""

# Import libraries
import os
import tarfile
import shutil 

import scanpy as sc
import anndata as ad

print(os.path.join(os.path.dirname(__file__), ".."))


project_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), ".."))

dir_endo    = os.path.join(project_dir, "data", "raw", "endo")
dir_cancer  =  os.path.join(project_dir, "data", "raw", "cancer")


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
    adatas, label="batch", keys=[p.rstrip("_") for p in prefixes]
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
]  # to remove doublets
adata_endo = adata_endo[
    adata_endo.obs.pct_counts_mt < 20, :
].copy()  # to remove cell with high mitochondrial activity

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
    adatas, label="batch", keys=[p.rstrip("_") for p in prefixes]
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
]  # to remove doublets
adata_cancer = adata_cancer[
    adata_cancer.obs.pct_counts_mt < 20, :
].copy()  # to remove cells with high mitochondrial activity

## Merge the two datasets
merged_adata = ad.concat(
    [adata_endo, adata_cancer],
    axis=0,  # Concatenate along cells
    join="inner",  # Keep intersected genes only
    keys=["Endometriosis", "Cancer"],
    merge="unique",  # Handle overlapping metadata uniquely
    index_unique="-",  # Avoid duplicate observation names
    fill_value=0,  # Fill missing values wth zeros
)



ad.AnnData.write_h5ad(merged_adata, filename="data/anndata.h5ad", compression='gzip')
