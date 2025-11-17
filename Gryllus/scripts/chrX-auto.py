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
    p = base / s / "outs"
    sp = p / "spatial"

    ad = sc.read_visium(path=p, count_file="filtered_feature_bc_matrix.h5")
    ad.obs["sample"] = s

    pos = pd.read_csv(sp / "tissue_positions.csv", header=None)
    pos.columns = ["barcode","in","row","col","px","py"]
    pos = pos.set_index("barcode").loc[ad.obs_names]

    ad.obs["px"] = pd.to_numeric(pos["py"], errors="coerce")
    ad.obs["py"] = pd.to_numeric(pos["px"], errors="coerce")

    img = Image.open(sp / "tissue_hires_image.png")
    sf = json.load(open(sp / "scalefactors_json.json"))["tissue_hires_scalef"]

    ad.obs["X_pixel"] = ad.obs["py"] * sf
    ad.obs["Y_pixel"] = ad.obs["px"] * sf
    ad.uns["img"] = np.array(img)

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
for s, ad in adatas.items():
    xp = [g for g in xg if g in ad.var_names]
    ap = [g for g in ag if g in ad.var_names]
    if not xp or not ap: continue

    Xs = ad[:, xp].X.sum(axis=1)
    As = ad[:, ap].X.sum(axis=1)
    if hasattr(Xs, "A1"):
        Xs = Xs.A1
        As = As.A1

    rat = Xs / (As + 1e-6)
    log = np.log1p(rat)

    ad.obs["X_ratio"] = rat
    ad.obs["log_X_ratio"] = log

    img = ad.uns["img"]

    # plot MAIN
    plt.figure(figsize=(6,6))
    plt.imshow(img)
    plt.scatter(ad.obs["X_pixel"], ad.obs["Y_pixel"], c=log, s=8, cmap="inferno")
    plt.axis("off")
    plt.show()

    # scatter
    plt.figure(figsize=(6,6))
    plt.scatter(auto_sum, X_sum, c=log_X_ratio, cmap="plasma", alpha=0.6, s=10)
    plt.xlabel("Autosome expression")
    plt.ylabel("chrX expression")
    plt.title(f"{s} – Scatter log-colored")
    plt.grid(True)
    plt.colorbar(fraction=0.046, pad=0.04, label="log(X/autosome)")
    plt.tight_layout()
    plt.show()

    # total
    plt.figure(figsize=(8,8))
    plt.imshow(img)
    plt.scatter(ad.obs["X_pixel"], ad.obs["Y_pixel"],
                c=ad.obs["total_expr"], cmap="viridis", s=10, alpha=0.8)
    plt.title(f"{s} – Total expression")
    plt.axis("off")
    plt.colorbar(fraction=0.046, pad=0.04, label="Total expression")
    plt.tight_layout()
    plt.show()

    # corr
    plt.figure(figsize=(6,6))
    plt.scatter(ad.obs["total_expr"], ad.obs["log_X_ratio"], alpha=0.5, s=10)
    plt.xlabel("Total expression")
    plt.ylabel("log(X/autosome)")
    plt.title(f"{s} – Corrélation total vs ratio")
    plt.grid(True)
    plt.tight_layout()
    plt.show()



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