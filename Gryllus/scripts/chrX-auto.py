import scanpy as sc
import pandas as pd
import numpy as np
import json
from pathlib import Path
from PIL import Image
import matplotlib.pyplot as plt
import gzip
import sys
sys.modules['tensorflow'] = None

# paths
base = Path("/home/martinpl/projects/datashare/grillus/IAB")
samples = ["SampleA_IAB", "SampleC_IAB", "SampleD_IAB"]

# load
adatas = {}
for s in samples:
    count_path = base / s / "outs"
    spatial_path = count_path / "spatial"

    ad = sc.read_visium(path=count_path, count_file="filtered_feature_bc_matrix.h5")
    ad.obs["sample"] = s

    positions = pd.read_csv(spatial_path / "tissue_positions.csv", header=None)
    positions.columns = ["barcode","in_tissue","array_row","array_col",
                         "pxl_row_in_fullres","pxl_col_in_fullres"]
    positions = positions.set_index("barcode").loc[ad.obs_names]

    ad.obs["pxl_row_in_fullres"] = pd.to_numeric(positions["pxl_row_in_fullres"], errors="coerce")
    ad.obs["pxl_col_in_fullres"] = pd.to_numeric(positions["pxl_col_in_fullres"], errors="coerce")

    img = Image.open(spatial_path / "tissue_hires_image.png")
    with open(spatial_path / "scalefactors_json.json") as f:
        scale_factor = json.load(f)["tissue_hires_scalef"]

    ad.obs["X_pixel"] = ad.obs["pxl_col_in_fullres"] * scale_factor
    ad.obs["Y_pixel"] = ad.obs["pxl_row_in_fullres"] * scale_factor

    ad.uns["img"] = np.array(img)      
    ad.uns["image_size"] = img.size

    adatas[s] = ad

# concat
adata = adatas["SampleA_IAB"].concatenate(
    adatas["SampleC_IAB"],
    adatas["SampleD_IAB"],
    batch_key="sample",
    batch_categories=samples
)

# map
gtf = base / "ref" / "genes.gtf.gz"
gene_chr = {}
with gzip.open(gtf, "rt") as f:
    for ln in f:
        if ln.startswith("#"): continue
        fs = ln.strip().split("\t")
        chrom = fs[0]
        info = fs[8]
        if 'gene_id "' not in info: continue
        gid = info.split('gene_id "')[1].split('"')[0]
        gene_chr.setdefault(gid, chrom)

genes = set(adata.var_names)
xg = [g for g in genes if gene_chr.get(g) == "chrX"]
ag = [g for g in genes if gene_chr.get(g,"").startswith("chr") and gene_chr.get(g)!="chrX"]

# log-ratio
# norm

for s, ad in adatas.items():
    xp = [g for g in xg if g in ad.var_names]
    ap = [g for g in ag if g in ad.var_names]
    if not xp or not ap: 
        continue

    # sums
    Xs = ad[:, xp].X.sum(axis=1)
    print(Xs)
    As = ad[:, ap].X.sum(axis=1)
    print(As)
    Ts = ad.X.sum(axis=1)

    # arr
    if hasattr(Xs, "A1"):
        Xs = Xs.A1
        As = As.A1
        Ts = Ts.A1

    # frac
    Xn = Xs / (Ts)
    An = As / (Ts)

    # ratio
    Rn = Xn / (An)
    #Ln = np.log1p(Rn)
    Ln = np.log(Rn)
    print(Ln)

    # save
    ad.obs["Xn"] = Xn
    ad.obs["An"] = An
    ad.obs["Xn_ratio"] = Rn
    ad.obs["log_Xn_ratio"] = Ln
    ad.obs["total_UMI"] = Ts

    img = ad.uns["img"]

    # plot MAIN
    plt.figure(figsize=(6,6))
    plt.imshow(img)
    sc = plt.scatter(ad.obs["X_pixel"], ad.obs["Y_pixel"], 
                    c=ad.obs["log_Xn_ratio"], s=8, cmap="inferno", alpha=0.9)
    plt.axis("off")

    cbar = plt.colorbar(sc, fraction=0.046, pad=0.04)
    cbar.set_label("log(Xn/An)", fontsize=9)

    plt.tight_layout()
    plt.show()


    # scatter
    plt.figure(figsize=(6,6))
    plt.scatter(An, Xn, c=Ln, cmap="plasma", alpha=0.6, s=10)
    plt.xlabel("Autosome (norm)")
    plt.ylabel("chrX (norm)")
    plt.title(f"{s} – Scatter log-normalized")
    plt.grid(True)
    plt.colorbar(fraction=0.046, pad=0.04, label="log(Xn/An)")
    plt.tight_layout()
    plt.show()

    # total
    plt.figure(figsize=(8,8))
    plt.imshow(img)
    plt.scatter(ad.obs["X_pixel"], ad.obs["Y_pixel"],
                c=ad.obs["total_UMI"], cmap="viridis", s=10, alpha=0.8)
    plt.title(f"{s} – Total UMI")
    plt.axis("off")
    plt.colorbar(fraction=0.046, pad=0.04, label="Total UMI")
    plt.tight_layout()
    plt.show()

    # corr
    plt.figure(figsize=(6,6))
    plt.scatter(ad.obs["total_UMI"], ad.obs["log_Xn_ratio"], alpha=0.5, s=10)
    plt.xlabel("Total UMI")
    plt.ylabel("log(Xn/An)")
    plt.title(f"{s} – Corr total vs log-ratio")
    plt.grid(True)
    plt.tight_layout()
    plt.show()


'''
# leiden clustering
for s, ad in adatas.items():
    try:
        # prep
        ad.raw = ad
        sc.pp.normalize_total(ad)
        sc.pp.log1p(ad)
        sc.pp.highly_variable_genes(ad, n_top_genes=1500, flavor="seurat", subset=True)
        sc.pp.scale(ad)
        sc.tl.pca(ad)
        sc.pp.neighbors(ad)
        sc.tl.leiden(ad, resolution=1.5, key_added="cl")

        # map
        xp = [g for g in xg if g in ad.var_names]
        ap = [g for g in ag if g in ad.var_names]

        out = []
        for c in ad.obs["cl"].unique():
            idx = ad.obs["cl"] == c
            sub = ad[idx,:]
            Xs = sub[:,xp].X.sum()
            As = sub[:,ap].X.sum()

            out.append({
                "cl": c,
                "X_ratio": Xs/(As+1e-6),
                "X_mean": sub[:,xp].X.mean(),
                "A_mean": sub[:,ap].X.mean()
            })

        df = pd.DataFrame(out)

        # barplot
        plt.figure(figsize=(8,3))
        plt.bar(df["cl"], df["X_ratio"])
        plt.tight_layout()
        plt.show()

        # heatmap
        h = df.set_index("cl")[["X_mean","A_mean"]]
        plt.figure(figsize=(6,3))
        sns.heatmap(h, annot=False)
        plt.tight_layout()
        plt.show()

    except Exception as e:
        print("err:", e)
        
'''
