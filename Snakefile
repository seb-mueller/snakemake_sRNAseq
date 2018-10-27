import pandas as pd
import os

# to make shell rule to work we need to determine the base path of Snakefile since we expect
# the scripts directory there as well
SRCDIR = srcdir("")

configfile: "config.yaml"
missmatches =  config['MAPPING']['missmatches']
reference   =  config['MAPPING']['reference']
refbase     =  os.path.basename(reference)
mode        =  config['MAPPING']['mode']
# threads     =  config['THREADS']

if mode == "unique":
    bowtie_par = "-m 1"
elif mode == "multi":
    bowtie_par = ""
else:
    bowtie_par = ""


# Get sample names from samples.csv
samples = pd.read_table("samples.csv", header=0, sep=',', index_col=0)

rule all:
    input:
        expand('logs/{sample}_R1_fastqc.html', sample=samples.index),
        # expand('trimmed/{sample}.fastq.gz', sample=samples.index),
        # expand('mapped/{sample}_{refbase}.bam', sample=samples.index, refbase=refbase),
        expand('mapped/{sample}.bam.bai', sample=samples.index, refbase=refbase),
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

rule bowtie:
    """maps small RNAs using bowtie and sorts them using samtools"""
    input:
        fastq="trimmed/{sample}.fastq.gz",
    output:
        "mapped/{sample}.bam"
    log:
        "logs/bowtie/{sample}.log"
    params:
        extra=""
    threads: 32
    shell:
        "bowtie {reference} --threads {threads} -v {missmatches} "
        "{bowtie_par} -q {input} -S 2> {log}"

rule postmapping:
    """bam.bai samtools flagstat etc"""
    input:
        "mapped/{sample}.bam"
    output:
        "mapped/{sample}.bam.bai"
    log:
        flagstat="logs/bowtie/{sample}_flagstat.log"
    shell:
        """
        samtools index {input}
        samtools flagstat {input} > {log.flagstat}
        """

 # bowtie [options]* <ebwt> {-1 <m1> -2 <m2> | --12 <r> |
 # --interleaved <i> | <s>} [<h
 #  -q                 query input files are FASTQ .fq/.fastq (default)
 #  -v <int>           report end-to-end hits w/ <=v mismatches; ignore qualities
 #  -p/--threads <int> number of alignment threads to launch (default: 1)
# Reporting:
 #  -k <int>           report up to <int> good alignments per read (default: 1)
 #  -a/--all           report all alignments per read (much slower than low -k)
 #  -m <int>           suppress all alignments if > <int> exist (def: no limit)
 #  -M <int>           like -m, but reports 1 random hit (MAPQ=0); requires --best
 #  --best             hits guaranteed best stratum; ties broken by quality
 #  --strata           hits in sub-optimal strata aren't reported (requires --best)
# samtools
  # -b       output BAM
  # -h       include header in SAM output

# bowtie data/index/genome --threads 1 -
# set -euo pipefail;  bowtie data/index/genome --threads 1 -v 0 -m 1 -q trimmed/test2.fastq.gz -S
# | samtools view -Sbh - > mapped/test2_genome.bam 2> logs/bowtie/test2_genome.log v 0 -m 1 -q trimmed/test3.fastq.gz | samtools view -Sbh - > mapped/test3_genome.bam
