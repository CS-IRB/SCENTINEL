import os
import scvelo as scv
import scanpy as sc
import numpy as np
import anndata
import loompy
import pandas as pd
import argparse

# Settings and directories
parser = argparse.ArgumentParser(description="Run scVelo analysis.")
parser.add_argument('--output_dir', required=True, help="Output directory for results.")
parser.add_argument('--merged_anndata_file', required=True, help="Path to the merged AnnData file.")
parser.add_argument('--loom_files', required=True, help="Comma-separated list of Loom file paths.")
parser.add_argument('--sample_names', required=True, nargs='+', help="List of sample names.")
args = parser.parse_args()

output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)

merged_anndata_file = args.merged_anndata_file

loom_files = {sample: path for sample, path in zip(args.sample_names, args.loom_files.split(','))}
print(f"Parsed loom_files: {loom_files}")

scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params('scvelo')

# Global lists to collect per-sample velocity matrices and AnnData objects.
combined_spliced, combined_unspliced = [], []
processed_anndatas = []

# Helper function to integrate velocity data for one sample
def add_velocity_to_anndata(merged_anndata, loom_file, sample_name):
    loom = loompy.connect(loom_file)
    
    # Get cell IDs from Loom (e.g. "sample:barcode") and clean them.
    loom_cells = loom.ca['CellID']
    loom_cells_cleaned = [cell.split(':')[1].rstrip('x') for cell in loom_cells]
    print(f"Loom cells (cleaned) for sample {sample_name}: {loom_cells_cleaned[:5]}... ({len(loom_cells_cleaned)} total)")
    
    # Get full cell names from AnnData for the given sample.
    merged_cells_full = merged_anndata.obs_names[merged_anndata.obs['sample'] == str(sample_name)].tolist()
    # Clean AnnData cell names (assuming format "barcode-...")
    merged_cells_cleaned = [cell.split('-')[0] for cell in merged_cells_full]
    print(f"Merged cells (cleaned) for sample {sample_name}: {merged_cells_cleaned[:5]}... ({len(merged_cells_cleaned)} total)")
    
    # Find common cleaned cell barcodes.
    common_cleaned = np.intersect1d(loom_cells_cleaned, merged_cells_cleaned)
    print(f"Common cells for sample {sample_name}: {len(common_cleaned)}")
    
    # Build list of full AnnData cell names that match the common cleaned barcodes.
    cells_to_keep = [cell_full for cell_full in merged_cells_full if cell_full.split('-')[0] in common_cleaned]
    
    if len(cells_to_keep) == 0:
        print(f"No matching cells found for sample {sample_name}.")
        loom.close()
        return None
    
    if len(merged_cells_full) > len(cells_to_keep):
        print(f"Subsetted merged AnnData from {len(merged_cells_full)} → {len(cells_to_keep)} cells for sample {sample_name}")
        merged_anndata = merged_anndata[cells_to_keep, :]
    
    # Subset Loom data to common cells.
    if len(loom_cells_cleaned) > len(common_cleaned):
        print(f"Subsetted Loom data from {len(loom_cells_cleaned)} → {len(common_cleaned)} cells for sample {sample_name}")
        loom_cells_idx = np.isin(loom_cells_cleaned, common_cleaned)
        loom_cells_idx_int = np.where(loom_cells_idx)[0]
        spliced = loom.layers['spliced'][:, loom_cells_idx_int]
        unspliced = loom.layers['unspliced'][:, loom_cells_idx_int]
    else:
        spliced = loom.layers['spliced']
        unspliced = loom.layers['unspliced']
    
    loom_genes = loom.ra['Gene']
    merged_genes = merged_anndata.var_names.tolist()
    # Clean loom gene names (e.g. remove suffixes)
    loom_genes_cleaned = [gene.split('_')[0] for gene in loom_genes]
    common_genes = np.intersect1d(loom_genes_cleaned, merged_genes)
    print(f"Common genes for sample {sample_name}: {len(common_genes)}")
    if len(common_genes) == 0:
        print(f"No common genes found for sample {sample_name}.")
        loom.close()
        return None
    
    # Subset Loom matrices to include only common genes.
    loom_genes_idx = np.isin(loom_genes_cleaned, common_genes)
    loom_genes_idx_int = np.where(loom_genes_idx)[0]
    spliced = spliced[loom_genes_idx_int, :]
    unspliced = unspliced[loom_genes_idx_int, :]
    
    # Subset AnnData to only include common genes.
    merged_anndata = merged_anndata[:, common_genes]
    
    # Transpose Loom matrices to shape (cells, genes)
    spliced = spliced.T
    unspliced = unspliced.T

    # Append matrices to global lists.
    combined_spliced.append(spliced)
    combined_unspliced.append(unspliced)
    
    loom.close()
    return merged_anndata


# Load the merged AnnData from file.
merged_anndata = anndata.read(merged_anndata_file)

# Process each sample separately.
processed_anndatas = []
for sample_name, loom_file in loom_files.items():
    subset_anndata = add_velocity_to_anndata(merged_anndata, loom_file, sample_name)
    if subset_anndata is not None:
        processed_anndatas.append(subset_anndata)
        print(f"Finished processing sample {sample_name} with {subset_anndata.n_obs} cells.")
    else:
        print(f"Skipping sample {sample_name} due to errors.")

# Concatenate the processed AnnDatas along the cell axis.
if len(processed_anndatas) == 0:
    raise ValueError("No samples were successfully processed!")
final_anndata = anndata.concat(processed_anndatas, join='outer', merge='same')

# Concatenate the corresponding spliced/unspliced matrices.
combined_spliced_mat = np.concatenate(combined_spliced, axis=0)
combined_unspliced_mat = np.concatenate(combined_unspliced, axis=0)
print(f"Final combined spliced shape: {combined_spliced_mat.shape}")
print(f"Final combined unspliced shape: {combined_unspliced_mat.shape}")
print(f"Final merged AnnData shape: {final_anndata.shape}")

# Validate that the final combined matrices match the final Anndata.
if combined_spliced_mat.shape[0] != final_anndata.n_obs:
    raise ValueError(f"Mismatch in cell count: spliced matrix has {combined_spliced_mat.shape[0]} cells but final AnnData has {final_anndata.n_obs}.")
if combined_spliced_mat.shape[1] != final_anndata.n_vars:
    raise ValueError("Mismatch in gene count between spliced matrix and final AnnData.")

# Assign the velocity layers.
final_anndata.layers['spliced'] = combined_spliced_mat
final_anndata.layers['unspliced'] = combined_unspliced_mat

# Save the updated merged AnnData.
final_anndata.write(os.path.join(output_dir, "merged_anndata_with_velocity.h5ad"))

# Downstream analysis functions
def perform_steady_state_analysis(adata, output_dir, sample_name=None):
    scv.pp.filter_and_normalize(adata, min_shared_counts=5, n_top_genes=3000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    scv.tl.velocity(adata)
    scv.tl.velocity_graph(adata)
    scv.tl.velocity_confidence(adata)
    scv.tl.velocity_pseudotime(adata)

    suffix = f"_{sample_name}" if sample_name else ""
    scv.pl.velocity_embedding_stream(adata, basis='umap', color='seurat_clusters',
                                     save=os.path.join(output_dir, f'velocity_embedding_stream{suffix}.pdf'))
    scv.tl.rank_velocity_genes(adata, groupby='seurat_clusters', min_corr=.3)
    pd.DataFrame(adata.uns['rank_velocity_genes']['names']).to_csv(
        os.path.join(output_dir, f'rank_velocity_genes{suffix}.csv')
    )
    scv.pl.velocity(adata, var_names=adata.var_names[:5], color='seurat_clusters',
                    save=os.path.join(output_dir, 'velocity_genes.pdf'))
    keys = 'velocity_length', 'velocity_confidence'
    scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95],
                   save=os.path.join(output_dir, 'scatter_velocity_length_confidence.pdf'))
    scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot',
                   save=os.path.join(output_dir, 'scatter_velocity_pseudotime.pdf'))
    adata.write(os.path.join(output_dir, f'scvelo_res_steady_state{suffix}.h5ad'))
    print(f"Finished steady-state velocity analysis for {'merged dataset' if not sample_name else sample_name}.")

def perform_dynamical_analysis(adata, output_dir, sample_name=None):
    scv.pp.filter_and_normalize(adata, min_shared_counts=5, n_top_genes=3000)
    sc.pp.pca(adata, n_comps=30)
    sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    scv.tl.recover_dynamics(adata)
    scv.tl.velocity(adata, mode='dynamical')
    scv.tl.velocity_graph(adata)
    scv.tl.velocity_pseudotime(adata)

    suffix = f"_{sample_name}" if sample_name else ""
    scv.pl.velocity_embedding_stream(adata, basis='umap', color='seurat_clusters',
                                     save=os.path.join(output_dir, f'dynamical_velocity_embedding_stream{suffix}.pdf'))
    scv.tl.rank_velocity_genes(adata, groupby='seurat_clusters', min_corr=.3)
    pd.DataFrame(adata.uns['rank_velocity_genes']['names']).to_csv(
        os.path.join(output_dir, f'dynamical_rank_velocity_genes{suffix}.csv')
    )
    scv.pl.velocity(adata, var_names=adata.var_names[:5], color='seurat_clusters',
                    save=os.path.join(output_dir, 'dynamical_velocity_genes.pdf'))
    keys = 'velocity_length', 'velocity_confidence'
    scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95],
                   save=os.path.join(output_dir, 'dynamical_scatter_velocity_length_confidence.pdf'))
    scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot',
                   save=os.path.join(output_dir, 'dynamical_scatter_velocity_pseudotime.pdf'))
    adata.write(os.path.join(output_dir, f'scvelo_res_dynamical{suffix}.h5ad'))
    print(f"Finished dynamical velocity analysis for {'merged dataset' if not sample_name else sample_name}.")

# Perform analyses on merged dataset
adata_merged = sc.read(os.path.join(output_dir, "merged_anndata_with_velocity.h5ad"), cache=True)
scv.pl.proportions(adata_merged, groupby='seurat_clusters',
                   save=os.path.join(output_dir, 'proportions.pdf'))
perform_steady_state_analysis(adata_merged, output_dir)
perform_dynamical_analysis(adata_merged, output_dir)

# Perform analyses on individual samples
for sample_name, loom_file in loom_files.items():
    print(f"Processing sample: {sample_name}")
    adata_sample = sc.read(os.path.join(output_dir, "merged_anndata_with_velocity.h5ad"), cache=True)
    adata_sample = adata_sample[adata_sample.obs['sample'] == str(sample_name)].copy()

    sample_output_dir = os.path.join(output_dir, f"sample_{sample_name}")
    os.makedirs(sample_output_dir, exist_ok=True)

    perform_steady_state_analysis(adata_sample, sample_output_dir, sample_name)
    perform_dynamical_analysis(adata_sample, sample_output_dir, sample_name)

print("Finished all analyses for merged dataset and individual samples.")
