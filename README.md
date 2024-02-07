
# Current Implementation


## Workflow Directory:

```bash
.
├── config
│   └── pipeline.yaml
├── envs
│   ├── cnv.yaml
│   ├── fusion.yaml
│   ├── metadata.yaml
│   ├── microarray.yaml
│   ├── mutation.yaml
│   ├── PharmacoSet.yaml
│   ├── rnaseq.yaml
│   ├── snakemake.yaml
│   ├── test.yaml
│   └── treatmentResponse.yaml
├── profiles
│   └── config.yaml
├── rules
│   ├── cnv.smk
│   ├── fusion.smk
│   ├── metadata.smk
│   ├── microarray.smk
│   ├── mutation.smk
│   ├── rnaseq.smk
│   └── treatmentResponse.smk
├── scripts
│   ├── build_PharmacoSet.R
│   ├── cnv
│   │   └── make_CNV_SE.R
│   ├── fusion
│   │   └── make_FUSION_SE.R
│   ├── metadata
│   │   ├── annotate_ChEMBL.R
│   │   ├── annotate_PubChemCIDS.R
│   │   ├── combine_annotated_treatmentData.R
│   │   ├── getCellosaurus
│   │   │   ├── getCellosaurusObject.R
│   │   │   └── mapCellosaurus.R
│   │   ├── map_treatments_to_PubChemCID.R
│   │   ├── preprocess_geneAnnotation.R
│   │   ├── preprocess_sampleMetadata.R
│   │   ├── preprocess_treatmentMetadata.R
│   │   └── utils.R
│   ├── microarray
│   │   └── make_MICROARRAY_SE.R
│   ├── mutation
│   │   └── make_MUTATION_SE.R
│   ├── rnaseq
│   │   └── make_RNASEQ_SE.R
│   ├── template.R
│   └── treatmentResponse
│       └── build_treatmentResponseExperiment.R
└── Snakefile
```  

# Setup

Ensure you have conda installed. If not, install [miniconda](https://docs.conda.io/en/latest/miniconda.html).
Install mamba for faster package management:

``` bash
conda install mamba -n base -c conda-forge 
```

Create a conda environment with the required dependencies:
> [!NOTE] 
> This pipeline uses `snakemake=7.32.4`. It is highly recommended to build the environment using the provided `snakemake.yaml` file.

```bash
mamba env create -f GDSC-Pharmacoset_Snakemake/workflow/envs/snakemake.yaml
```

Activate the environment:

```bash
conda activate gdsc_pharmacoset
```


## Workflow Execution:

### Main Run
> [!NOTE]
> The following command was tested on a linux machine, with 30 cores and 128GB of RAM.
> The pipeline can be run on a machine with less resources, but it will take longer to complete.
> To modify the number of cores used, change the `--cores` flag to the desired number.

``` bash
snakemake --snakefile workflow/Snakefile \
    --cores 30 \
    --show-failed-logs \
    --keep-going \
    --rerun-incomplete \
    --printshellcmds \
    --use-conda \
    --use-singularity
```


``` bash
snakemake --profile workflow/profiles/
```

### Dry Run
```bash
snakemake --profile workflow/profiles/ --dryrun
```

### Create all Conda Environments
> [!TIP] 
> Creating each conda environment can take a long time when running the entire pipeline. 
> To create all conda environments without running the pipeline, use the `--conda-create-envs-only` flag.

```bash
snakemake --profile workflow/profiles/ --use-conda --conda-create-envs-only
```


## RuleGraph 
``` bash
snakemake --profile workflow/profiles/ --rulegraph | dot -Tsvg > resources/rulegraph.svg
```

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

