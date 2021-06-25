import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")


##### load config and sample sheets #####

configfile: "config.yaml"

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
aligners=config["params"]["aligners"].split(",")
print("Aligners:", aligners)


##### target rules #####

rule all:
    input:
        "qc/multiqc_report_pretrim.html",
        "qc/multiqc_report_posttrim.html",
        expand("qc/multiqc_report_{aligner}.html", aligner=aligners),
        expand("results/{aligner}/featureCounts/all.fixed.featureCounts", aligner=aligners)


include: "rules/qc.smk"
include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/featurecounts.smk"
