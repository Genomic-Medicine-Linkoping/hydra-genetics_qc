__author__ = "Martin R"
__copyright__ = "Copyright 2021, Martin R"
__email__ = "martin.rippin@scilifelab.uu.se"
__license__ = "GPL-3"


rule mosdepth:
    input:
        bam="alignment/merge_bam/{sample}_{type}.bam",
        bai="alignment/merge_bam/{sample}_{type}.bam.bai",
    output:
        bed="qc/mosdepth/{sample}_{type}.regions.bed.gz",
        csi="qc/mosdepth/{sample}_{type}.regions.bed.gz.csi",
        glob="qc/mosdepth/{sample}_{type}.mosdepth.global.dist.txt",
        region="qc/mosdepth/{sample}_{type}.mosdepth.region.dist.txt",
        summary="qc/mosdepth/{sample}_{type}.mosdepth.summary.txt",
    params:
        by=config.get("mosdepth", {}).get("by", ""),
        extra=config.get("mosdepth", {}).get("extra", ""),
        bed=config.get("mosdepth", {}).get("bed", ""),
    log:
        "qc/mosdepth/{sample}_{type}.log",
    benchmark:
        repeat(
            "qc/mosdepth/{sample}_{type}.benchmark.tsv",
            config.get("mosdepth", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("mosdepth", {}).get("threads", config["default_resources"]["threads"])
    resources:
        threads=config.get("mosdepth", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("mosdepth", {}).get("time", config["default_resources"]["time"]),
        mem_mb=config.get("mosdepth", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("mosdepth", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("mosdepth", {}).get("partition", config["default_resources"]["partition"]),
    container:
        config.get("mosdepth", {}).get("container", config["default_container"])
    conda:
        "../envs/mosdepth.yaml"
    message:
        "{rule}: Calculating coverage for {wildcards.sample}_{wildcards.type}"
    wrapper:
        "0.80.2/bio/mosdepth"
