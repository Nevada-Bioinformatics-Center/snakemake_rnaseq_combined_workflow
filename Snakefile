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
    sample="[\w\\-\\.]+",
    unit="rep\d+",
    trimmer="[a-z]+"

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
aligners=config["params"]["aligners"].split(",")
trimmers=config["params"]["trimmers"].split(",")
pese=config["params"]["pese"]
print("Aligners:", aligners)
print("Trimmers:", trimmers)
print("PE/SE mode:", pese)
cwd = os.getcwd()
print("Cwd:", cwd)

def strip_suffix(pattern, suffix):
    return pattern[: -len(suffix)]


wrappers_version="v1.9.0"

##### target rules #####
rule all:
    input:
        expand("qc/multiqc_report_pretrim_{pese}.html", pese=pese),
        #expand("qc/multiqc_report_posttrim_{trimmer}.html", trimmer=trimmers),
        expand("qc/multiqc_report_{aligner}_{trimmer}_{pese}.html", aligner=aligners, trimmer=trimmers, pese=pese),
        #expand("qc/multiqc_report_{aligner}_{trimmer}_nofct.html", aligner=aligners, trimmer=trimmers),
        expand(cwd+"/results/{aligner}/all.{aligner}.{trimmer}_{pese}.fixcol2.featureCounts", aligner=aligners, trimmer=trimmers, pese=pese),
        expand(cwd+"/results/{aligner}/all.{aligner}.{trimmer}_{pese}_multi.fixcol2.featureCounts", aligner=aligners, trimmer=trimmers, pese=pese),
        expand(cwd+"/results/{aligner}/all.{aligner}.{trimmer}_{pese}_multifrac.fixcol2.featureCounts", aligner=aligners, trimmer=trimmers, pese=pese),


include: "rules/qc.smk"
include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/featurecounts.smk"
