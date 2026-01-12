import os
import pandas as pd
import numpy as np
import scanpy as sc
import SpaGCN as spg
import random, torch
import warnings
import matplotlib.pyplot as plt
import cv2
import anndata

warnings.filterwarnings("ignore")

# Base paths
base_visium = "/home/martinpl/projects/datashare/PDAC/visium_PDAC"
out_dir     = "/home/martinpl/projects/PDAC/spagcn_output"

os.makedirs(out_dir, exist_ok=True)

slides = [
    "Visium_FFPE_V43T08-051_D", "Visium_FFPE_V44L22-378_D", "Visium_FFPE_V44L23-391_D",
    "Visium_FFPE_V43T08-041_A", "Visium_FFPE_V44L01-325_A", "Visium_FFPE_V44L23-362_A",
    "Visium_FFPE_V43T08-041_D", "Visium_FFPE_V44L01-325_D", "Visium_FFPE_V44L23-362_D",
    "Visium_FFPE_V43T08-051_A", "Visium_FFPE_V44L22-378_A", "Visium_FFPE_V44L23-391_A"
]

# Helper: find spatial folder (outs/outs_old)
def get_spatial_input(input_dir):
    for subdir in ["outs", "outs_old"]:
        spatial_path = os.path.join(input_dir, subdir)
        if os.path.isdir(spatial_path):
            return spatial_path
    raise FileNotFoundError(f"{input_dir} does not contain 'outs' or 'outs_old'")


for sample_name in slides:
    print(f"\n--- Processing {sample_name} ---")

    sample_dir   = os.path.join(base_visium, sample_name)
    spatial_input = get_spatial_input(sample_dir)
    filtered_data_input = os.path.join(sample_dir, f"{sample_name}_filtered.csv")

    # Output dir per slide
    out_slide = os.path.join(out_dir, sample_name)
    os.makedirs(out_slide, exist_ok=True)

    # CSV -> h5ad 
    df = pd.read_csv(filtered_data_input, index_col=0)
    adata = anndata.AnnData(X=df.values)
    adata.obs_names = df.index
    adata.var_names = df.columns
    adata.write_h5ad(f"{out_slide}/{sample_name}_filtered_normalized.h5ad")

    # Load spatial info 
    spatial = pd.read_csv(
        f"{spatial_input}/spatial/tissue_positions.csv",
        sep=",", header=0, na_filter=False, index_col=0
    )

    adata.obs["x1"] = spatial["in_tissue"]
    adata.obs["x2"] = spatial["array_row"]
    adata.obs["x3"] = spatial["array_col"]
    adata.obs["x4"] = spatial["pxl_row_in_fullres"]
    adata.obs["x5"] = spatial["pxl_col_in_fullres"]

    adata.obs["x_array"] = adata.obs["x2"]
    adata.obs["y_array"] = adata.obs["x3"]
    adata.obs["x_pixel"] = adata.obs["x4"]
    adata.obs["y_pixel"] = adata.obs["x5"]

    adata.var_names = [i.upper() for i in list(adata.var_names)]
    adata.var["genename"] = adata.var.index.astype("str")
    adata.write_h5ad(f"{out_slide}/{sample_name}_sample_data.h5ad")
    adata = sc.read(f"{out_slide}/{sample_name}_sample_data.h5ad")

    # Histology image 
    img = cv2.imread(f"{spatial_input}/spatial/cytassist_image.tiff")

    x_pixel = adata.obs["x_pixel"].astype(float).round().astype(int).values
    y_pixel = adata.obs["y_pixel"].astype(float).round().astype(int).values

    # Adjacency 
    s = 1
    b = 49
    adj = spg.calculate_adj_matrix(
        x=x_pixel, y=y_pixel,
        x_pixel=x_pixel, y_pixel=y_pixel,
        image=img, beta=b, alpha=s, histology=True
    )
    np.savetxt(f"{out_slide}/{sample_name}_adj.csv", adj, delimiter=",")

    adata = sc.read(f"{out_slide}/{sample_name}_sample_data.h5ad")
    adj = np.loadtxt(f"{out_slide}/{sample_name}_adj.csv", delimiter=",")
    adata.var_names_make_unique()

    # SpaGCN parameters 
    p = 0.5
    l = spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)

    r_seed = t_seed = n_seed = 100
    random.seed(r_seed)
    torch.manual_seed(t_seed)
    np.random.seed(n_seed)

    res = 0.5
    clf = spg.SpaGCN()
    clf.set_l(l)

    clf.train(adata, adj, init_spa=True, init="louvain",
              res=res, tol=5e-3, lr=0.05, max_epochs=200)
    y_pred, prob = clf.predict()
    adata.obs["pred"] = y_pred
    adata.obs["pred"] = adata.obs["pred"].astype('category')
    adata.write_h5ad(f"{out_slide}/{sample_name}_results.h5ad")

    # Plot 
    plot_color = [
        "#F56867","#FEB915","#C798EE","#59BE86","#7495D3","#D1D1D1",
        "#6D1A9C","#15821E","#3A84E6","#997273","#787878","#DB4C6C",
        "#9E7A7A","#554236","#AF5F3C","#93796C","#F9BD3F","#DAB370",
        "#877F6C","#268785"
    ]

    domains = "pred"
    num_celltype = len(adata.obs[domains].unique())
    adata.uns[domains+"_colors"] = list(plot_color[:num_celltype])

    ax = sc.pl.scatter(
        adata, alpha=1, x="y_pixel", y="x_pixel",
        color=domains, title=domains, show=False,
        size=100000/adata.shape[0]
    )
    ax.set_aspect("equal", "box")
    ax.axes.invert_yaxis()
    plt.savefig(f"{out_slide}/{sample_name}_pred.png", dpi=600)
    plt.close()

    print(f"✓ {sample_name} finished -> {out_slide}")