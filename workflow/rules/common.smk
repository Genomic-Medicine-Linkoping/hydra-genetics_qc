# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas A"
__copyright__ = "Copyright 2021, Jonas A"
__email__ = "jonas.almlof@igp.uu.se"
__license__ = "GPL-3"

import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *

min_version("6.8.0")

### Set and validate config file


configfile: "config.yaml"


validate(config, schema="../schemas/config.schema.yaml")
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")


### Read and validate samples file

samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

### Read and validate units file

units = (
    pandas.read_table(config["units"], dtype=str)
    .set_index(["sample", "type", "flowcell", "lane", "barcode"], drop=False)
    .sort_index()
)
validate(units, schema="../schemas/units.schema.yaml")


def get_flowcell(units, wildcards):
    flowcells = set([u.flowcell for u in get_units(units, wildcards)])
    if len(flowcells) > 1:
        raise ValueError("Sample type combination from different sequence flowcells")
    return flowcells.pop()


### Set wildcard constraints


wildcard_constraints:
    sample="|".join(samples.index),
    unit="N|T|R",


def compile_output_list(wildcards):
    output_files = [
        "qc/hotspot_info/%s_%s.hotspot_info.tsv" % (sample, unit_type)
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
    ]
    output_files.append("qc/multiqc/multiqc.html")
    return output_files


def get_unit_for_fastqc(units: pandas.DataFrame, wildcards: snakemake.io.Wildcards) -> pandas.Series:
    """
    function used to extract one unit(row) from units.tsv
    Args:
        units: DataFrame generate by importing a file following schema defintion
               found in pre-alignment/workflow/schemas/units.schema.tsv
        wildcards: wildcards object with at least the following wildcard names
               sample, type, flowcell, and lane
    Returns:
        Series containing data of the selected row
    Raises:
        raises an exception (KeyError) if no unit can be extracted from the Dataframe
    """
    unit = units.loc[(wildcards.sample, wildcards.type)].dropna()
    return unit


def get_fastq_file_for_fastqc(units: pandas.DataFrame, wildcards: snakemake.io.Wildcards, read_pair: str = "fastq1") -> str:
    """
    function used to extract path for one unit(row) from units.tsv
    Args:
        units: DataFrame generate by importing a file following schema defintion
               found in pre-alignment/workflow/schemas/units.schema.tsv
        wildcards: wildcards object with at least the following wildcard names
               sample, type, flowcell, and lane
        read_pair: fast1 or fastq2
    Returns:
        path for fastq file as a str
    Raises:
        raises an exception (KeyError) if no unit can be extracted from the Dataframe
    """
    unit = get_unit_for_fastqc(units, wildcards)
    if read_pair not in ["fastq1", "fastq2"]:
        raise ValueError("Incorrect input value error {}: expected {} or {}".format(read_pair, "fastq1", "fastq2"))
    return unit[read_pair]
