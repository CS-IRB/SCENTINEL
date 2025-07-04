import os
import subprocess
import argparse
import pandas as pd

# Parse arguments
parser = argparse.ArgumentParser(description="Run pySCENIC pipeline.")
parser.add_argument('--output_dir', required=True, help="Output directory")
parser.add_argument('--loom_file', required=True, help="Input loom file with expression data")
parser.add_argument('--tf_list', required=True, help="List of transcription factors")
parser.add_argument('--db_glob', required=True, help="Path to cisTarget databases")
parser.add_argument('--n_tasks', type=int, default=4, help="Number of parallel tasks")
parser.add_argument('--annotations', required=True, help="Path to annotations file")
args = parser.parse_args()

# Define file paths
loom_file = args.loom_file
tf_list = args.tf_list
db_glob = args.db_glob
output_dir = args.output_dir
n_tasks = args.n_tasks
annotations = args.annotations

# Step 1: Gene Regulatory Network (GRN) Inference
grn_output = os.path.join(output_dir, 'adjacencies.tsv')
subprocess.run([
    'pyscenic', 'grn', '--method', 'grnboost2',
    '--num_workers', str(n_tasks),
    '--output', grn_output,
    loom_file,
    tf_list
], check=True)

# Step 2: Regulon Prediction (cisTarget Analysis)
regulons_output = os.path.join(output_dir, 'regulons.csv')
motifs_output = os.path.join(output_dir, 'motifs.csv')
subprocess.run([
    'pyscenic', 'ctx',
    '--annotations_fname', annotations,
    '--expression_mtx_fname', loom_file,
    '--output', regulons_output,
    '--mask_dropouts',
    '--num_workers', str(n_tasks),
    grn_output,
    db_glob
], check=True)

# Step 3: Cellular Enrichment (AUCell Analysis)
auc_mtx_output = os.path.join(output_dir, 'auc_mtx.csv')
loom_output = os.path.join(output_dir, 'output.loom')
subprocess.run([
    'pyscenic', 'aucell',
    '--output', loom_output,
    '--num_workers', str(n_tasks),
    loom_file,
    regulons_output
], check=True)

# def create_sif_file(grn_output, regulons_output, output_dir, 
#                    interaction_type="Tf-target", 
#                    importance_threshold=3.0, 
#                    min_regulon_size=0,
#                    expressed_tfs=None):
#     """Create filtered SIF file enforcing both importance and regulon filters"""
#     try:
#         # 1. Load and filter adjacencies
#         adjacencies = pd.read_csv(grn_output, sep='\t', header=0)
#         adjacencies.columns = ['TF', 'Target', 'Importance']
        
#         print(f"Original adjacencies: {len(adjacencies):,} edges")
#         adjacencies = adjacencies[adjacencies['Importance'] > importance_threshold]
#         print(f"After importance filter (> {importance_threshold}): {len(adjacencies):,} edges")
        
#         if expressed_tfs is not None:
#             adjacencies = adjacencies[adjacencies['TF'].isin(expressed_tfs)]
#             print(f"After TF expression filter: {len(adjacencies):,} edges")

#         # 2. Load regulons and calculate target counts
#         regulons_df = pd.read_csv(regulons_output)
#         regulons_df['nTargets'] = regulons_df['TargetGenes'].apply(
#             lambda x: len(str(x).split(',')) if pd.notna(x) else 0
#         )
#         regulons_df['nTargets'] = regulons_df['nTargets'] // 2

#         print(f"\nOriginal regulons: {len(regulons_df):,}")
#         regulons_df = regulons_df[regulons_df['nTargets'] >= min_regulon_size]
#         print(f"After size filter (≥ {min_regulon_size} targets): {len(regulons_df):,} regulons")

#         # 3. Create set of valid TF-target pairs from filtered adjacencies
#         valid_pairs = set(zip(adjacencies['TF'], adjacencies['Target']))
#         print(f"\nValid TF-target pairs passing all filters: {len(valid_pairs):,}")

#         # 4. Strict filtering - only keep edges that exist in both sources
#         edges = set()
#         invalid_count = 0
        
#         for _, row in regulons_df.iterrows():
#             tf = row['TF']
#             targets_str = str(row['TargetGenes'])
            
#             try:
#                 clean_str = (targets_str.replace('frozenset({', '')
#                               .replace('})', '')
#                               .replace('[', '')
#                               .replace(']', '')
#                               .replace('(', '')
#                               .replace(')', '')
#                               .replace("'", ""))
                
#                 for item in clean_str.split(','):
#                     item = item.strip()
#                     if not item or item.replace('.', '').isdigit():
#                         continue
                        
#                     target = ''.join(c for c in item if c.isalnum() or c in ['-', '_'])
                    
#                     if len(target) > 1:
#                         if (tf, target) in valid_pairs:  # KEY FILTER
#                             edges.add((tf, target))
#                         else:
#                             invalid_count += 1
#                     else:
#                         invalid_count += 1
                        
#             except Exception as e:
#                 print(f"Skipping {tf}: {e}")
#                 continue

#         # 5. Save final filtered SIF
#         sif_path = os.path.join(output_dir, 'strict_filtered_network.sif')
#         with open(sif_path, 'w') as f:
#             f.write("\n".join([f"{tf}\t{interaction_type}\t{target}" for tf, target in edges]))
        
#         # 6. Final report
#         print(f"\n=== FINAL NETWORK ===")
#         print(f"Edges passing all filters: {len(edges):,}")
#         print(f"TFs in network: {len({tf for tf,_ in edges}):,}")
#         print(f"Target genes: {len({tg for _,tg in edges}):,}")
#         print(f"Invalid entries skipped: {invalid_count:,}")
#         print(f"\nSaved strict filtered SIF to: {sif_path}")
        
#         return sif_path
        
#     except Exception as e:
#         print(f"Error: {e}")
#         raise


# create_sif_file(grn_output, regulons_output, output_dir)