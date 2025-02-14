
# Installing and running Cell2Location + cluster execution

## Option 1: Conda environment 

### Setting up the Conda Environment on the Cluster

To build the required Conda environment, use the `environment.yml` configuration file available in the official repository: 

`wget https://github.com/BayraktarLab/cell2location/blob/master/environment.yml`

Given the resource-intensive nature of environment creation, it is recommended to submit a job to the cluster for building the environment:

```bash
#!/bin/bash

source /applis/environments/conda.sh

conda env create --prefix /home/user/cell2loc_env -f /user/cell2loc_env.yml

```

### Resolving Dependency Conflicts

After creating the environment, additional manual adjustments may be required to ensure compatibility with the latest version of Cell2Location and associated packages:

Uptdated `.yml` file: `(git link to be updated):***`

Manually update specific packages within the environment to avoid compatibility issues:

```bash
conda activate cell2loc_env
pip install anndata==0.10.8
pip install scvi-tools==1.0.4
pip install jaxlib==0.4.23
pip install jax==0.4.23
pip install numpy==1.26.4
pip install scipy==1.12.0

```
### Aditionnal step por A100 GPU (bigfoot default)

```bash
export PYTHONNOUSERSITE="True"
conda create -y -n cell2location_cuda118_torch22 python=3.10
conda activate cell2location_cuda118_torch22

pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

pip3 install scvi-tools==1.1.2

pip install git+https://github.com/BayraktarLab/cell2location.git#egg=cell2location[tutorials,dev]
python -m ipykernel install --user --name=cell2location_cuda118_torch22 --display-name='Environment (cell2location_cuda118_torch22)'
```



## Option 2: Singularity (TO UPDATE)

To utilize Cell2Location via Singularity, download the container image `.sif` of the latest version:

```bash
wget https://cell2location.cog.sanger.ac.uk/singularity/cell2location-v0.06-alpha.sif
```

Submit your computational jobs to the OAR scheduler using the `singularity exec` command.

Note: The latest available container version (`cell2location-v0.06-alpha.sif`) contains outdated dependencies that could cause compatibility issues. Refer to the following issue for more details:

https://github.com/BayraktarLab/cell2location/issues/137

# Running Cell2Location

For guidance on running Cell2Location, refer to the official tutorial: https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html#Loading-packages

Python script for testing: `(git link to be updated):***`

1. Update the paths in your script to point to your results directory and set a symbolic link to `/bettik`

```python
results_folder = '/home/User/cell2location/results/'
```

2. Modify the `batch_size` variable regarding your GPU power (default settings will lead to memory usage problems, `batch_size=500` gave a correct execution of the test)


3. To load packages for cell2location execution I had to use this modification of the orignal code:

```python
import sys
import subprocess
import os
IN_COLAB = "google.colab" in sys.modules
if IN_COLAB:
    subprocess.run(["pip", "install", "--quiet", "scvi-colab"], check=True)
    from scvi_colab import install
    install()
    subprocess.run(["pip", "install", "--quiet", "git+https://github.com/BayraktarLab/cell2location#egg=cell2location[tutorials]"], check=True)
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cell2location
from matplotlib import rcParams
```

#### Execution in cluster Dahu (CPU)

While CPU execution is possible, it is generally discouraged due to the performance limitations. If you opt to use a CPU, ensure that you modify the `accelerator_connector.py` to use the CPU backend to avoid compatibility issues. A tutorial can be followed here: https://github.com/Lightning-AI/pytorch-lightning/issues/13991

Submit OAR job via:

```bash
oarsub --project (project name) -p "network_address='(luke node)'" -l nodes=1,walltime=48:00:00 /home/user/cell2location/run_test.sh

```

#### Execution in cluster Bigfoot (GPU)

Due to the high memory demands and image processing of Cell2Location, it is strongly recommended to use a GPU for execution. Note that the Bigfoot cluster may have long queue times, so it is advisable to debug your Python script using a minimal core allocation before scaling up to full GPU usage.

Remember to create a symbolic link to your `/bettik` repo in the results folder in order to avoid storage problems in the cluster.

```bash
oarsub -l /nodes=1/core=1/gpu=1 -p "gpumodel='A100'" --project (project name)) -S ./run_test_env.sh
```

This configuration should allow Cell2Location to complete its execution successfully:

```bash
oarsub -l /nodes=1/core=32/gpu=1,walltime=48:00:00 -p "gpumodel='A100' and mem_per_gpu > 64" --project (project name) -S ./run_test.sh

```
