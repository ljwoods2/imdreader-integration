Full walkthrough
================

Create IMDReader mamba env
--------------------------

Assuming you're in a login node to start, get a compute node:
```bash
salloc -p general
```

From the compute node, create env with name "imdreader-test":
```bash
mkdir -p workspace
cd workspace
git clone https://github.com/Becksteinlab/imdreader.git
cd imdreader
module load mamba/latest
mamba env create --file devtools/conda-envs/test_env.yaml
```

Run the script
--------------
Clone repo and run script:
```bash
cd ~/workspace
git clone https://github.com/ljwoods2/imdreader-integration.git
cd slurm_scripts/gmx_imdreader
```

