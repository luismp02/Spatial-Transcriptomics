```markdown
# Tutorial: Installing `cell2location` with GPU support on ARM64 (GH200 Grace Hopper bigfoot Node)
This guide walks you through setting up a complete environment for spatial transcriptomics using `cell2location`, compiled from source to support GPU execution on an **ARM64 (aarch64)** Grace Hopper (GH200) architecture with an Nvidia H100 GPU.
Tested on:
- Cluster: `bigfoot.ciment`
- Node: `bigfoot-gh1`
- OS: Ubuntu aarch64
- CUDA: 12.1 (preinstalled)
- PyTorch: built from source, version `nightly`

## Step 1 — Launch interactive job on Grace Hopper node and activate conda
```bash
oarsub -I -t gh -l /nodes=1,walltime=2:00:00 --project=pr-epimed
```
If you haven’t installed Miniconda for aarch64 yet, [get it here](https://docs.conda.io/en/latest/miniconda.html).
Once installed:
```bash
source ~/miniconda3/bin/activate
```
Create a clean environment for building PyTorch:
```bash
conda create -n cell2loc_arm python=3.10 -y
conda activate cell2loc_arm
```

## Step 2 — Install build dependencies and clone Pytorch (from conda, no sudo needed)
```bash
conda install -c conda-forge cmake=3.26 ninja libprotobuf protobuf
conda install -c conda-forge setuptools cffi typing_extensions future six requests
```
> Note: although we installed `ninja`, we ended up using `make` due to compatibility issues with `libstdc++`.

```bash
git clone --recursive https://github.com/pytorch/pytorch
cd pytorch
git checkout nightly
git pull
git submodule sync
git submodule update --init --recursive
```

## Step 3 — Configure environment for CUDA build
```bash
export USE_CUDA=1
export USE_CUDNN=0
export USE_CUSPARSELT=0
export BUILD_TEST=0
export TORCH_CUDA_ARCH_LIST="8.9"
export CMAKE_GENERATOR="Unix Makefiles"
export CC=/usr/bin/gcc
export CXX=/usr/bin/g++
export PATH=/usr/local/cuda/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH
nvcc --version  # should work
```
## Step 4 — Build PyTorch from source
```bash
python setup.py clean
rm -rf build
python setup.py bdist_wheel
```
```bash
cd dist
pip install torch-*.whl
cd ~
```
Test installation:
```bash
python -c "import torch; print(torch.__version__); print(torch.cuda.is_available())"
# → Should print the torch version and True
```

## Step 5 — Install `cell2location` with GPU support
```bash
pip install "cell2location[GPU]"
pip install scanpy anndata pandas matplotlib seaborn jax tqdm h5py
pip uninstall torchaudio torchvision  # optional cleanup
```
- Use `/var/tmp/$USER` for fast NVMe scratch storage.
```
