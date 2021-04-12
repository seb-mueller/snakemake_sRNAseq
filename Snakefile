import pandas as pd
import os

# Usage:
# source activate srna_mapping
# snakemake --use-conda --conda-prefix ~/.myconda -n

# to make shell rule to work we need to determine the base path of Snakefile since we expect
# the scripts directory there as well
SRCDIR = srcdir("")

configfile: "config.yaml"
missmatches =  config['MAPPING']['missmatches']
reference   =  config['MAPPING']['reference']
refbase     =  os.path.basename(reference)
# mode        =  config['MAPPING']['mode']
multi_m     =  config['MAPPING']['m']
extra_params_bowtie = config['MAPPING']['extra_params_bowtie'] 
# threads     =  config['THREADS']

if multi_m == 1:
  mode = "unique"
elif multi_m > 1:
  mode = "multi" + str(multi_m)
else:
    raise Exception('m should be specified (>=1) in config.yaml. The value of m was: {}'.format(x))


# Get sample names from samples.csv
samples = pd.read_table("samples.csv", header=0, sep=',', index_col=0)

rule all:
    input:
        expand('logs/fastqc/raw/{sample}_fastqc.html', sample=samples.index),
        expand('logs/fastqc/trimmed/{sample}_fastqc.html', sample=samples.index),
        # expand('trimmed/{sample}.fastq.gz', sample=samples.index),
        expand('mapped/{sample}_MappedOn_{refbase}_{mode}.bam.bai', sample=samples.index, refbase=refbase, mode=mode),
        # expand('mapped/{sample}.bam.bai', sample=samples.index, refbase=refbase),
        # 'reports/fastqc.html',
        expand('mapped/bws/{sample}_MappedOn_{refbase}_{mode}.cpm.bw', sample=samples.index, refbase=refbase, mode=mode),

rule trimall:
    input:
        expand('logs/fastqc/raw/{sample}_fastqc.html', sample=samples.index),
        expand('logs/fastqc/trimmed/{sample}_fastqc.html', sample=samples.index),
        expand('trimmed/{sample}.fastq.gz', sample=samples.index),

rule fq2fa:
    input:
        expand('trimmed/fasta/{sample}.fa', sample=samples.index),

rule staridx:
    input:
      #"staridx"
      config['STAR']['star_idx_dir']
      #expand('{reference}_star', reference=reference)

rule starmap:
    input:
        expand('mapped_star/{sample}.Aligned.sortedByCoord.out.bam', sample=samples.index)

rule fastqc_raw:
    """Create fastqc report"""
    input:
        'data/{sample}.fastq.gz'
    output:
        html='logs/fastqc/raw/{sample}_fastqc.html',
        zip= 'logs/fastqc/raw/{sample}_fastqc.zip'
    params: '--extract'
    log:
        "logs/fastqc/raw/{sample}.log"
    wrapper:
        '0.70.0/bio/fastqc'

# rule multiqc_fastqc_reads:
#     """https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/multiqc.html"""
#     input:
#         expand('logs/{sample}_fastqc.html', sample=samples.index)
#     output:
#         html='reports/fastqc.html',
#         txt='reports/fastqc_data/multiqc_general_stats.txt'
#     params: '-m fastqc'
#     wrapper:
#         '0.49.0/bio/multiqc'

if config['FILTER']['cutadapt']['quality-filter'] != 0:
  qual = " -q " + str(config['FILTER']['cutadapt']['quality-filter'])
else:
  qual = ""

rule cutadapt:
    input:
        'data/{sample}.fastq.gz'
    output:
        fastq="trimmed/{sample}.fastq.gz",
        qc="trimmed/{sample}.qc.txt"
    params:
      adapters=" -a " +     config['FILTER']['cutadapt']['adapter'],
        extra=
          " -m " + str(config['FILTER']['cutadapt']['minimum-length']) +
          " -M " + str(config['FILTER']['cutadapt']['maximum-length']) +
          " "    + str(config['FILTER']['cutadapt']['extra-params']) +
          qual
    log:
        "logs/cutadapt/{sample}.log"
    threads: 1 # set desired number of threads here
    wrapper:
        "0.70.0/bio/cutadapt/se"

rule fastq_to_fasta:
    """Convert fq 2 fa"""
    input:
        'trimmed/{sample}.fastq.gz'
    output:
        'trimmed/fasta/{sample}.fa',
    shell:
        """
        fastq_to_fasta -i <(zcat {input}) -o {output}
        """

rule fastqc_trimmed:
    """Create fastqc report"""
    input:
        'trimmed/{sample}.fastq.gz'
    output:
        html='logs/fastqc/trimmed/{sample}_fastqc.html',
        zip= 'logs/fastqc/trimmed/{sample}_fastqc.zip'
    params: '--extract'
    log:
        "logs/fastqc/trimmed/{sample}.log"
    wrapper:
        '0.70.0/bio/fastqc'

rule star_index:
    input:
        fasta = config['MAPPING']['reference']
    output:
        #directory("{reference}_star")
        directory(config['STAR']['star_idx_dir'])
    message:
        "Creating STAR index in for " + 
        config['MAPPING']['reference'] + " in " + 
        config['STAR']['star_idx_dir']
    threads:
        1
    params:
        extra = ""
    log:
        "logs/star_index.log"
        # "logs/star_index_{reference}.log"
    wrapper:
        "0.70.0/bio/star/index"

rule star_se:
    input:
        fq1 = "trimmed/{sample}.fastq.gz",
    output:
        # see STAR manual for additional output files
        # "star/{sample}_{refbase}/Aligned.out.sam"
        "mapped_star/{sample}.Aligned.sortedByCoord.out.bam",
    log:
        "logs/star/{sample}.log"
    params:
        # path to STAR reference genome index
        # index="staridx",
        index=config['STAR']['star_idx_dir'],
    threads: 8
    shell:
      "STAR "
      #"{extra} "
      "--runThreadN {threads} "
      "--genomeDir {params.index} "
      "--readFilesIn {input.fq1} "
      "--readFilesCommand zcat "
      " --runMode "               + config["STAR"]["runMode"] +
      " --genomeLoad "            + config["STAR"]["genomeLoad"] +
      " --outSAMtype "            + config["STAR"]["outSAMtype"] +
      " --alignEndsType "         + config["STAR"]["alignEndsType"] +
      " --outFileNamePrefix "     + "mapped_star/{wildcards.sample}." +
      " --scoreDelOpen "          + str(config["STAR"]["scoreDelOpen"]) +
      " --scoreInsOpen "          + str(config["STAR"]["scoreInsOpen"]) +
      " --alignIntronMax "        + str(config["STAR"]["alignIntronMax"]) +
      " --outFilterMismatchNmax " + str(config["STAR"]["outFilterMismatchNmax"]) +
      " --outFilterMultimapNmax " + str(config["STAR"]["outFilterMultimapNmax"]) +
      " --alignSJDBoverhangMin "  + str(config["STAR"]["alignSJDBoverhangMin"]) +
      " --limitBAMsortRAM "       + str(config["STAR"]["limitBAMsortRAM"]) +
      "--outStd Log "
      "{log}"  
      #"0.70.0/bio/star/align" # didn't work since outFileNamePrefix couldn't be overwritten

rule bowtie:
    """maps small RNAs using bowtie and sorts them using samtools"""
    input:
        fastq="trimmed/{sample}.fastq.gz",
    output:
        "mapped/{sample}_MappedOn_{refbase}_{mode}.bam"
    log:
        "logs/bowtie/{sample}_MappedOn_{refbase}_{mode}.log"
    params:
        # this is adding read-group ID to bam header, e.g.:
        # @RG     ID:SRR1234_leaf
        # SRR10419098.5387676     16      chr01  ... RG:Z:SRR10419098_leaf NH:i:7  HI:i:1  n
        # can also be acchieve with picard tools
        # picard AddOrReplaceReadGroups I=$file O=for_shortstack/rg_${file} RGID=${file%%.bam} \

        "--sam-RG ID:{sample}"
    threads: 16
    conda: 'environment.yaml'
    shell:
        "bowtie {reference} --threads {threads} -v {missmatches}"
        " -m {multi_m} {extra_params_bowtie} {params} -q {input} -S 2> {log}"
        " | samtools view -Sbh - | samtools sort -o {output}"

rule postmapping:
    """bam.bai samtools flagstat etc"""
    input:
        "mapped/{sample}.bam"
    output:
        "mapped/{sample}.bam.bai"
    log:
        flagstat     = "logs/bowtie/{sample}_flagstat.log",
        idxstat      = "logs/bowtie/{sample}_idxstats.log",
        coveragestat = "logs/bowtie/{sample}_coveragestats.log",
        depthstat    = "logs/bowtie/{sample}_depthstats.log"
    conda: 'environment.yaml'
    shell:
        """
        samtools index        {input}
        samtools flagstat     {input} >    {log.flagstat}
        samtools idxstats     {input} >    {log.idxstat}
        # bamtools coverage -in {input} -out {log.coveragestat}
        # above gives:Pileup::Run() : Data not sorted correctly!
        # bamtools depth -in    {input} -out {log.depthstat}
        """

rule calccoverage:
    """compute coverage using deeptools into bigWig(bw) file"""
    input:
        bam="mapped/{sample}.bam",
        bai="mapped/{sample}.bam.bai"
    output:
        "mapped/bws/{sample}.cpm.bw"
    params:
        binsize="10"
    threads: 8
    conda: 'environment.yaml'
    shell:
        """
        bamCoverage -b {input.bam} -o {output}  --normalizeUsing CPM --binSize {params.binsize} -p {threads}
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
