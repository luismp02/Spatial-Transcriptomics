#!/usr/bin/env python

"""
Author : @mgautheret
"""

import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
from scipy.sparse import issparse
import random, torch
import warnings
warnings.filterwarnings("ignore")
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import SpaGCN as spg
import cv2
import sys
import anndata

## == Set parameters == ##
sample_name="Visium_FFPE_V43T08-051_D"
sample_dir="/home/genouest/cnrs_umr6290/mgauth/real_data/"
out_dir="/home/genouest/cnrs_umr6290/mgauth/out/outSpaGCN/"



## == Convert csv to h5ad and save h5ad file == ##
df = pd.read_csv(f"{sample_dir}{sample_name}/filtered_data/{sample_name}_filtered_normalized.csv", index_col=0)
adata = anndata.AnnData(X=df.values)
adata.obs_names = df.index
adata.var_names = df.columns
adata.write_h5ad(f'{out_dir}{sample_name}_filtered_normalized.h5ad')



## == Load data == ##
spatial=pd.read_csv(f"{sample_dir}{sample_name}/outs/spatial/tissue_positions.csv",sep=",",header=0,na_filter=False,index_col=0) 

adata.obs["x1"] = spatial["in_tissue"]
adata.obs["x2"] = spatial["array_row"]
adata.obs["x3"] = spatial["array_col"]
adata.obs["x4"] = spatial["pxl_row_in_fullres"]
adata.obs["x5"] = spatial["pxl_col_in_fullres"]

adata.obs["x_array"]=adata.obs["x2"]
adata.obs["y_array"]=adata.obs["x3"]
adata.obs["x_pixel"]=adata.obs["x4"]
adata.obs["y_pixel"]=adata.obs["x5"]

#Select captured samples
adata.var_names=[i.upper() for i in list(adata.var_names)]
adata.var["genename"]=adata.var.index.astype("str")
adata.write_h5ad(f"{out_dir}{sample_name}_sample_data.h5ad")
adata=sc.read(f"{out_dir}{sample_name}_sample_data.h5ad")

#Read in hitology image
img=cv2.imread(f"{sample_dir}{sample_name}/outs/spatial/cytassist_image.tiff")

#Set coordinates
x_array=adata.obs["x_array"].tolist()
y_array=adata.obs["y_array"].tolist()
x_pixel=adata.obs["x_pixel"].tolist()
y_pixel=adata.obs["y_pixel"].tolist()

x_pixel = adata.obs["x_pixel"].astype(float).round().astype(int).values
y_pixel = adata.obs["y_pixel"].astype(float).round().astype(int).values

## == Calculate adjacency matrix == ##
s=1 #The s parameter determines the weight given to histology when calculating Euclidean distance between every two spots. ‘s = 1’ means that the histology pixel intensity value has the same scale variance as the (x,y) coordinates, whereas higher value of ‘s’ indicates higher scale variance, hence, higher weight to histology, when calculating the Euclidean distance. 
b=49 # The b parameter determines the area of each spot when extracting color intensity.

adj=spg.calculate_adj_matrix(x=x_pixel,y=y_pixel, x_pixel=x_pixel, y_pixel=y_pixel, image=img, beta=b, alpha=s, histology=True)
#If histlogy image is not available, SpaGCN can calculate the adjacent matrix using the function below
#adj=spg.alculate_adj_matrix(x=x_pixel,y=y_pixel, histology=False)
np.savetxt(f'{out_dir}{sample_name}_adj.csv', adj, delimiter=',')


adata=sc.read(f"{out_dir}{sample_name}_sample_data.h5ad")
adj=np.loadtxt(f'{out_dir}{sample_name}_adj.csv', delimiter=',')
adata.var_names_make_unique()


p=0.5 #Percentage of total expression contributed by neighborhoods.
l=spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100) #Parameter to control p.


#Set seed
r_seed=t_seed=n_seed=100
random.seed(r_seed)
torch.manual_seed(t_seed)
np.random.seed(n_seed)

#If the number of clusters known, we can use the spg.search_res() fnction to search for suitable resolution(optional)
#n_clusters=5

#Search for suitable resolution, if n_clusters is defined
#res=spg.search_res(adata, adj, l, n_clusters, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)
res=0.5
clf=spg.SpaGCN()
clf.set_l(l)



## == Run clustering == ##
clf.train(adata,adj,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)
y_pred, prob=clf.predict()
adata.obs["pred"]= y_pred
adata.obs["pred"]=adata.obs["pred"].astype('category')

adata.write_h5ad(f"{out_dir}{sample_name}_results.h5ad")

adata=sc.read(f"{out_dir}{sample_name}_results.h5ad")


## == Plot clustering results == ##
#Set colors used
plot_color=["#F56867","#FEB915","#C798EE","#59BE86","#7495D3","#D1D1D1","#6D1A9C","#15821E","#3A84E6","#997273","#787878","#DB4C6C","#9E7A7A","#554236","#AF5F3C","#93796C","#F9BD3F","#DAB370","#877F6C","#268785"]

#Plot spatial domains
domains="pred"
num_celltype=len(adata.obs[domains].unique())
adata.uns[domains+"_colors"]=list(plot_color[:num_celltype])
ax=sc.pl.scatter(adata,alpha=1,x="y_pixel",y="x_pixel",color=domains,title=domains,color_map=plot_color,show=False,size=100000/adata.shape[0])
ax.set_aspect('equal', 'box')
ax.axes.invert_yaxis()
plt.savefig(f"{out_dir}{sample_name}_pred.png", dpi=600)
plt.close()


