# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule fastqc:
    input:
        "prealignment/merged/{sample}_{type}_{read}.fastq.gz",
    output:
        html=temp("qc/fastqc/{sample}_{type}_{read}_fastqc.html"),
        zip=temp("qc/fastqc/{sample}_{type}_{read}_fastqc.zip"),
        tmp=temp("qc/fastqc/{sample}_{type}_{read}/tmp.txt"),
    params:
        outdir=lambda wildcards, output: os.path.dirname(os.path.abspath(output.html)),
        tmp=lambda wildcards, output: os.path.dirname(os.path.abspath(output.tmp)),
    log:
        "qc/fastqc/{sample}_{type}_{read}.fastqc.log",
    benchmark:
        repeat("qc/fastqc/{sample}_{type}_{read}.fastqc.benchmark.tsv", config.get("fastqc", {}).get("benchmark_repeats", 1))
    threads: config.get("fastqc", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("fastqc", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("fastqc", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("fastqc", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("fastqc", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("fastqc", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("fastqc", {}).get("container", config["default_container"])
    conda:
        "../envs/fastqc.yaml"
    message:
        "{rule}: qc/{rule}/{wildcards.sample}_{wildcards.type}_{wildcards.read}"
    shell:
        "(touch {output.tmp} && "
        "fastqc "
        "--quiet "
        "-t {threads} "
        "-d {params.tmp} "
        "--outdir {params.outdir} "
        "{input}) "
        "&> {log}"
