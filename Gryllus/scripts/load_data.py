import argparse
import anndata
import pandas as pd
import scipy.io
from scipy.io import mmread

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outs", required=True, help="Space Ranger outs directory (contains filtered_feature_bc_matrix and spatial)")
    ap.add_argument("--files", required=True, help="Directory with signature CSVs (liste_gene_and_code.csv etc.)")
    ap.add_argument("--out_h5ad", required=True, help="Output .h5ad path")
    args = ap.parse_args()

    matrix_path  = f"{args.outs}/filtered_feature_bc_matrix"
    spatial_path = f"{args.outs}/spatial"
    data_path    = args.files

    matrix = scipy.io.mmread(f"{matrix_path}/matrix.mtx.gz").tocsc()
    barcodes = pd.read_csv(f"{matrix_path}/barcodes.tsv.gz", header=None, names=["barcode"])
    features = pd.read_csv(f"{matrix_path}/features.tsv.gz", header=None, names=["gene_id", "gene_name", "feature_type"])

    genes_growth    = pd.read_csv(f"{data_path}/de_gene_names_growth.csv")
    genes_regen_up  = pd.read_csv(f"{data_path}/ur_regeneration.csv")
    genes_regen_down= pd.read_csv(f"{data_path}/dr_regeneration.csv")
    genes_scarring  = pd.read_csv(f"{data_path}/gene_name_scaring.csv")

    features['gene_id'] = (features['gene_id']
        .astype(str)
        .str.replace(r'\t', '', regex=True)
        .str.replace('Gene Expression', '', regex=True)
        .str.strip().str.lower()
        .str.slice(0, 9)
    )

    for df in [genes_growth, genes_regen_up, genes_regen_down, genes_scarring]:
        df[df.columns[0]] = (df[df.columns[0]].astype(str)
            .str.replace(r'\t', '', regex=True)
            .str.strip().str.lower()
        )

    adata = anndata.AnnData(X=matrix.T, obs=barcodes, var=features.set_index('gene_id'))
    adata.var.index = adata.var.index.str.slice(0, 9)

    # Nota: tu código tenía un typo: growth_gene usaba genes_scarring
    adata.var['growth_gene']     = adata.var.index.isin(genes_growth[genes_growth.columns[0]].tolist())
    adata.var['regen_up_gene']   = adata.var.index.isin(genes_regen_up[genes_regen_up.columns[0]].tolist())
    adata.var['regen_down_gene'] = adata.var.index.isin(genes_regen_down[genes_regen_down.columns[0]].tolist())
    adata.var['scarring_gene']   = adata.var.index.isin(genes_scarring[genes_scarring.columns[0]].tolist())

    tissue_positions = pd.read_csv(f"{spatial_path}/tissue_positions.csv")
    adata.obs = adata.obs.merge(tissue_positions, on="barcode", how="left")

    adata.write_h5ad(args.out_h5ad)

if __name__ == "__main__":
    main()