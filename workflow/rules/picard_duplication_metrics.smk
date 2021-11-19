# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2021, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule picard_duplication_metrics:
    input:
        bam="alignment/merge_bam/{sample}_{type}.bam",
    output:
        metrics=temp("qc/picard_duplication_metrics/{sample}_{type}.duplication_metrics.txt"),
    log:
        "qc/picard_duplication_metrics/{sample}_{type}.log",
    benchmark:
        repeat(
            "qc/picard_duplication_metrics/{sample}_{type}.benchmark.tsv",
            config.get("picard_duplication_metrics", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("picard_duplication_metrics", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("picard_duplication_metrics", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("picard_duplication_metrics", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("picard_duplication_metrics", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("picard_duplication_metrics", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("picard_duplication_metrics", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("picard_duplication_metrics", {}).get("container", config["default_container"])
    conda:
        "../envs/picard_duplication_metrics.yaml"
    message:
        "{rule}: Calculate qc using picard: qc/{rule}/{wildcards.sample}_{wildcards.type}"
    shell:
        "(picard CollectDuplicateMetrics INPUT={input.bam} M={output.metrics}) &> {log}"
