import pandas as pd

# to make shell rule to work we need to determine the base path of Snakefile since we expect
# the scripts directory there as well
SRCDIR = srcdir("")

configfile: "config.yaml"

# Get sample names from samples.csv
samples = pd.read_table("samples.csv", header=0, sep=',', index_col=0)

# a pseudo-rule that collects the target files
# rule all:
#     input:
#          expand(config["OUTDIR"] + "/{sample}_R1.fastq.gz", sample=config["SAMPLES"])

rule all:
    input:
        expand('logs/{sample}_R1_fastqc.html', sample=samples.index),
        expand('trimmed/{sample}.fastq.gz', sample=samples.index),
        'reports/fastqc.html',

rule fastqc:
    """Create fastqc report"""
    input:
        'data/{sample}_R1.fastq.gz'
    output:
        html='logs/{sample}_R1_fastqc.html',
        zip='logs/{sample}_R1_fastqc.zip'
    params: '--extract'
    log:
        "logs/fastqc/{sample}.log"
    wrapper:
        '0.27.1/bio/fastqc'

rule multiqc_fastqc_reads:
    """https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/multiqc.html"""
    input:
        expand('logs/{sample}_R1_fastqc.html', sample=samples.index)
    output:
        html='reports/fastqc.html',
        txt='reports/fastqc_data/multiqc_general_stats.txt'
    params: '-m fastqc'
    wrapper:
        '0.27.1/bio/multiqc'

rule cutadapt:
    input:
        'data/{sample}_R1.fastq.gz'
    output:
        fastq="trimmed/{sample}.fastq.gz",
        qc="trimmed/{sample}.qc.txt"
    params:
        " -a " +     config['FILTER']['cutadapt']['adapter'] +
        " -q " + str(config['FILTER']['cutadapt']['quality-filter']) +
        " -m " + str(config['FILTER']['cutadapt']['minimum-length']) +
        " -M " + str(config['FILTER']['cutadapt']['maximum-length'])
    log:
        "logs/cutadapt/{sample}.log"
    wrapper:
        "0.27.1/bio/cutadapt/se"
