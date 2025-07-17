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
    sample=r"[\w\.-]+",
    unit=r"rep\d+",
    trimmer=r"[a-z]+"

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


wrappers_version="v7.0.0"

# Build input file lists conditionally before rule
multiqc_pretrim = expand("qc/multiqc_report_pretrim_{pese}.html", pese=pese)

aligner_inputs = []
if "salmon" not in aligners:
    aligner_inputs = expand(
        cwd + "/results/{aligner}/all.{aligner}.{trimmer}_{pese}.fixcol2.featureCounts",
        aligner=aligners, trimmer=trimmers, pese=pese
    )

salmon_inputs = []
if "salmon" in aligners:
    salmon_inputs = expand(
        "salmon/{trimmer}_{pese}/{unit.sample}.{unit.unit}/quant.sf",
        trimmer=trimmers, pese=pese, unit=units.itertuples()
    )

##### target rules #####
rule all:
    input:
        multiqc_pretrim,
        expand("qc/multiqc_report_{aligner}_{trimmer}_{pese}.html", aligner=aligners, trimmer=trimmers, pese=pese),
        aligner_inputs,
        salmon_inputs
#        expand("qc/multiqc_report_pretrim_{pese}.html", pese=pese),
#        (
#            expand("qc/multiqc_report_{aligner}_{trimmer}_{pese}.html", aligner=aligners, trimmer=trimmers, pese=pese),
#            expand(cwd+"/results/{aligner}/all.{aligner}.{trimmer}_{pese}.fixcol2.featureCounts", aligner=aligners, trimmer=trimmers, pese=pese),
#            if "salmon" not in aligners else []
#        ),
#        (
#            expand("salmon/{trimmer}_{pese}/{unit.sample}.{unit.unit}/quant.sf", trimmer=trimmers, pese=pese, unit=units.itertuples())
#            if "salmon" in aligners else []
#        ),
#        #expand("qc/multiqc_report_{aligner}_{trimmer}_nofct.html", aligner=aligners, trimmer=trimmers),
#        #expand(cwd+"/results/{aligner}/all.{aligner}.{trimmer}_{pese}_multi.fixcol2.featureCounts", aligner=aligners, trimmer=trimmers, pese=pese),
#        #expand(cwd+"/results/{aligner}/all.{aligner}.{trimmer}_{pese}_multifrac.fixcol2.featureCounts", aligner=aligners, trimmer=trimmers, pese=pese),


include: "rules/qc.smk"
include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/featurecounts.smk"
