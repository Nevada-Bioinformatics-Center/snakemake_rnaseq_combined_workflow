import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")


##### load config and sample sheets #####

configfile: "config.yaml"

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
aligners=config["params"]["aligners"].split(",")
trimmers=config["params"]["trimmers"].split(",")
#samples=units["sample"]
print("Aligners:", aligners)
print("Trimmers:", trimmers)
#print("Samples:", samples)

def strip_suffix(pattern, suffix):
    return pattern[: -len(suffix)]


##### target rules #####
rule all:
    input:
        "qc/multiqc_report_pretrim.html",
        expand("qc/multiqc_report_posttrim_{trimmer}.html", trimmer=trimmers),
        expand("qc/multiqc_report_{aligner}_{trimmer}.html", aligner=aligners, trimmer=trimmers),
        #expand("results/{aligner}/featureCounts/all.fixed.featureCounts", aligner=aligners)


include: "rules/qc.smk"
include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/featurecounts.smk"