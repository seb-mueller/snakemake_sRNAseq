# snakemake_sRNAseq
Snakemake workflow for processing small RNA-seq libaries

# creating conda environment
```
conda env create --file environment.yaml --name srna_mapping
```

# activate 

```
source activate srna_mapping
```

# Usage:

```
snakemake --use-conda --conda-prefix ~/.myconda
```
