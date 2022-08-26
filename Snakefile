import pandas as pd
from snakemake.utils import validate, min_version
import os
##### set minimum snakemake version #####
min_version("5.1.2")


##### load config and sample sheets #####


#WRAPPER_PREFIX='https://github.com/hans-vg/snakemake-wrappers/raw'
WRAPPER_PREFIX='https://raw.githubusercontent.com/hans-vg/snakemake-wrappers'

configfile: "config.yaml"

wildcard_constraints:
    sample="[\w-]+",
    trimmer="\w+"

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
aligners=config["params"]["aligners"].split(",")
trimmers=config["params"]["trimmers"].split(",")
print("Aligners:", aligners)
print("Trimmers:", trimmers)
cwd = os.getcwd()
print("Cwd:", cwd)

def strip_suffix(pattern, suffix):
    return pattern[: -len(suffix)]

wrappers_version="v1.9.0"

##### target rules #####
rule all:
    input:
        "qc/multiqc_report_pretrim.html",
        #expand("qc/multiqc_report_posttrim_{trimmer}.html", trimmer=trimmers),
        expand("qc/multiqc_report_{aligner}_{trimmer}.html", aligner=aligners, trimmer=trimmers),
        #expand(cwd+"/results/{aligner}/all.{aligner}.{trimmer}.fixcol2.featureCounts", aligner=aligners, trimmer=trimmers),


include: "rules/qc.smk"
include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/featurecounts.smk"
