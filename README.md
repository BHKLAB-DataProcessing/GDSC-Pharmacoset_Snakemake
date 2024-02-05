
# Current Implementation

## Directed Acyclic Graph (DAG)
``` 
snakemake --profile workflow/profiles/ --dag | dot -Tsvg > resources/dag.svg
```

![DAG](./resources/dag.svg)

## filegraph
```
snakemake --profile workflow/profiles/  --filegraph | dot -Tsvg > resources/filegraph.svg
```

![filegraph](./resources/filegraph.svg)


# Setup

Ensure you have conda installed. If not, install [miniconda](https://docs.conda.io/en/latest/miniconda.html).
Install mamba for faster package management:

``` bash
conda install mamba -n base -c conda-forge 
```

Create a conda environment with the required dependencies:

```bash
mamba env create -f GDSC-Pharmacoset_Snakemake/workflow/envs/snakemake.yaml
```

Activate the environment:

```bash
conda activate gdsc_pharmacoset
```