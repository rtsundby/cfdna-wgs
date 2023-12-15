##################################################################
###   Integration testing snakefile for WGS cfDNA Processing   ###
##################################################################

import pandas as pd
import re
import numpy as np

configfile: "config/nf1_hg19.yaml"


# Values directly from configuration file
threads = config["threads"]
cfdna_wgs_threads = config["threads"]
default_container = config["default_container"]
cfdna_wgs_container = config["cfdna_wgs_container"]
genome_fasta = config["genome_fasta"]
genome_ref = config["genome_ref"]
cfdna_wgs_repo = config["cfdna_wgs_repo"]
adapters = config["adapters"]

# Directory values derived from datadir in configuration YAML
datadir                   = config["datadir"]
outdir                    = config["outdir"]
DOWNSAMPLE                = config["downsample"]
cfdna_wgs                 = datadir + "/analysis/cfdna_wgs"
cfdna_wgs_bams            = outdir + "/cfdna-wgs-bams"
cfdna_wgs_fastqs          = datadir + "/analysis/cfdna_wgs/fastqs"
qcdir                     = datadir + "/analysis/qc"
benchdir                  = datadir + "/benchmark"
logdir                    = datadir + "/logs"
refdir                    = config["ref_dir"] 
rawfastq                  = config["raw_fastq"]         

cfdna_wgs_scriptdir = config["cfdna_wgs_repo"] +  "/scripts"

###   Functions   ###
# Setup sample name index as a python dictionary
cfdna_wgs_libraries = pd.read_table(config["samples"])
print(cfdna_wgs_libraries)
# Make the dictionary
cfdna_wgs_library_indict = cfdna_wgs_libraries["library"].tolist()
cfdna_wgs_file_indict = cfdna_wgs_libraries["file"].tolist()
cfdna_wgs_lib_dict = dict(zip(cfdna_wgs_library_indict, cfdna_wgs_file_indict))

for k,v in cfdna_wgs_lib_dict.items():
    print(k,v)

CFDNA_WGS_LIBRARIES = list(cfdna_wgs_lib_dict.keys())
CFDNA_WGS_FASTQS = list(cfdna_wgs_lib_dict.values())

### smk workflow ###

rule all:
    input:
         expand(cfdna_wgs_fastqs +
            "/{library}_{processing}_{read}.fastq.gz",
             library = CFDNA_WGS_LIBRARIES,
             processing = ["processed", "unpaired"],
             read = ["R1", "R2"]),
         expand(cfdna_wgs_bams + "/{library}_{processing}.bam",
             library = CFDNA_WGS_LIBRARIES,
             processing = ["raw", "dedup", "filt"]),
         expand(cfdna_wgs_bams + "/{library}_ds{downsample}.bam",
             library = CFDNA_WGS_LIBRARIES,
             downsample = DOWNSAMPLE),

rule cfdna_wgs_fastp:
    benchmark: benchdir + "/{library}_cfdna_wgs_fastp.benchmark.txt",
    resources:
        time   = "8:00:00",
        mem_gb = "15g",
        cpus   = "10",
    input:
        read1 = cfdna_wgs_fastqs + "/{library}_raw_R1.fastq.gz",
        read2 = cfdna_wgs_fastqs + "/{library}_raw_R2.fastq.gz",
    log:
        cmd = logdir + "/{library}_cfdna_wgs_fastp.log",
        html = logdir + "/{library}_cfdna_wgs_fastp.html",
        json = logdir + "/{library}_cfdna_wgs_fastp.json",
    output:
        read1 = outdir + "/{library}_processed_R1.fastq.gz",
        read2 = outdir + "/{library}_processed_R2.fastq.gz",
        failed = outdir + "/{library}_failed_fastp.fastq.gz",
        unpaired1 = outdir + "/{library}_unpaired_R1.fastq.gz",
        unpaired2 = outdir + "/{library}_unpaired_R2.fastq.gz",
    params:
        script = cfdna_wgs_scriptdir + "/fastp_plos_pkg.sh",
        threads = cfdna_wgs_threads,
    shell:
        """
        {params.script} \
        {input.read1} \
        {input.read2} \
        {log.html} \
        {log.json} \
        {output.read1} \
        {output.read2} \
        {output.failed} \
        {output.unpaired1} \
        {output.unpaired2} \
        {resources.cpus} &> {log.cmd}
        """

# Align reads with BWA
rule cfdna_wgs_align:
    benchmark: benchdir + "/{library}_cfdna_wgs_align.benchmark.txt",
    resources:
        time   = "10:00:00",
        mem_gb = "100g",
        cpus   = "50",
    input:
        ref = genome_ref,
        read1 = outdir + "/{library}_processed_R1.fastq.gz",
        read2 = outdir + "/{library}_processed_R2.fastq.gz",
    log: logdir + "/{library}_cfdna_wgs_align.log",
    output:
        sort = cfdna_wgs_bams + "/{library}_raw.bam",
        index = cfdna_wgs_bams + "/{library}_raw.bam.bai",
    params:
        script = cfdna_wgs_scriptdir + "/align_plos_pkg.sh",
    shell:
        """
        {params.script} \
        {input.ref} \
        {input.read1} \
        {input.read2} \
        {resources.cpus} \
        {output.sort} &> {log}
        """

# Remove PCR duplicates from aligned reads
rule cfdna_wgs_dedup:
    benchmark: benchdir + "/{library}_cfdna_wgs_dedup.benchmark.txt",
    resources:
        time   = "8:00:00",
        mem_gb = "80g",
        cpus   = "8",
    input: cfdna_wgs_bams + "/{library}_raw.bam",
    log: logdir + "/{library}_cfdna_wgs_dedup.log",
    output: cfdna_wgs_bams + "/{library}_dedup.bam",
    params:
        script = cfdna_wgs_scriptdir + "/dedup_plos_pkg.sh",
        threads = cfdna_wgs_threads,
    shell:
        """
        {params.script} \
        {input} \
        {output} \
        {resources.cpus} &> {log}
        """

# Filter de-duplicated alignments.
# Remove unmapped, not primary, and duplicate reads. Additional location filter by config bedfile variable.

checkpoint cfdna_wgs_filter_alignment:
    benchmark: benchdir + "/{library}_cfdna_wgs_filter_alignment.benchmark.txt",
    resources:
        time   = "9:00:00",
        mem_gb = "2g",
        cpus   = "10",
    input: cfdna_wgs_bams + "/{library}_dedup.bam",
    log: logdir + "/{library}_cfdna_wgs_filter_alignment.log",
    output: cfdna_wgs_bams + "/{library}_filt.bam",
    params:
        script = cfdna_wgs_scriptdir + "/filter_alignment_plos_pkg.sh",
        threads = cfdna_wgs_threads,
    shell:
        """
        {params.script} \
        {input} \
        {resources.cpus} \
        {output} &> {log}
        """

rule downsample_bams:
    resources:
        time   = "9:00:00",
        mem_gb = "80g",
        cpus   = "10"
    input: cfdna_wgs_bams + "/{library}_filt.bam",
    #output: touch(logdir + "/{library}_{downsample}_downsample.done"),
    output: cfdna_wgs_bams + "/{library}_ds{downsample}.bam"
    params:
        out_dir = cfdna_wgs_bams,
        script = cfdna_wgs_scriptdir + "/downsample_bams_plos_pkg.sh",
        threads = cfdna_wgs_threads,
        suffix = "_filt.bam"
    shell:
        """
        {params.script} \
        {input} \
        {wildcards.downsample} \
        {params.out_dir} \
        {params.suffix} \
        {params.threads}
        """
#Get read quality by FASTQC
rule cfdna_wgs_fastqc:
    benchmark: benchdir+ "/{library}_{processing}_{read}_cfdna_wgs_fastqc.benchmark.txt",
    resources:
        time   = "2:00:00",
        mem_gb = "1g",
        cpus   = "2",
    input: cfdna_wgs_fastqs + "/{library}_{processing}_{read}.fastq.gz",
    log: logdir + "/{library}_{processing}_{read}_cfdna_wgs_fastqc.log",
    output:
        qcdir + "/{library}_{processing}_{read}_fastqc.html",
        qcdir + "/{library}_{processing}_{read}_fastqc.zip",
    params:
        outdir = qcdir,
        script = cfdna_wgs_scriptdir + "/fastqc.sh",
        threads = cfdna_wgs_threads,
    shell:
        """
        {params.script} \
        {input} \
        {params.outdir} \
        {resources.cpus} &> {log}
        """

# Get alignment QC using samtools
rule cfdna_wgs_alignment_qc:
    resources:
        time   = "2:00:00",
        mem_gb = "1g",
        cpus   = "2",
    input: cfdna_wgs_bams + "/{library}_{processing}.bam",
    log:
        flagstat = logdir + "/{library}_{processing}_flagstat_cfdna_wgs_alignment_qc.log",
        samstat = logdir + "/{library}_{processing}_samstats_cfdna_wgs_alignment_qc.log",
    output:
        flagstat = qcdir + "/{library}_{processing}_flagstat.txt",
        samstat = qcdir + "/{library}_{processing}_samstats.txt",
    params:
        script = cfdna_wgs_scriptdir + "/alignment_qc.sh",
        threads = cfdna_wgs_threads,
    shell:
        """
        {params.script} \
        {input} \
        {log.flagstat} \
        {log.samstat} \
        {output.flagstat} \
        {output.samstat} \
        {resources.cpus}
        """