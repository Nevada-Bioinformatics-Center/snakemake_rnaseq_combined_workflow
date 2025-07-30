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
    pese=r"pe|se",
    trimmer=r"[a-z]+"

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
has_fq2 = units["fq2"].notna() & (units["fq2"] != "")

if has_fq2.all():
    pese = "pe"
elif (~has_fq2).all():
    pese = "se"
else:
    raise ValueError(
        "Mixed single-end and paired-end samples detected in units.tsv.  "
        "Either every row needs BOTH fq1+fq2 (for PE) or NONE should have fq2 (for SE)."
    )


aligners=config["params"]["aligners"].split(",")
trimmers=config["params"]["trimmers"].split(",")
runmulti=config["params"]["runfctmulti"]
runmultifrac=config["params"]["runfctmultifrac"]
print("Aligners:", aligners)
print("Trimmers:", trimmers)
print("PE/SE mode:", pese)
cwd = os.getcwd()
print("Cwd:", cwd)

def strip_suffix(pattern, suffix):
    return pattern[: -len(suffix)]


wrappers_version="v7.0.0"

# 1) partition aligners
non_salmon = [a for a in aligners if a.lower() != "salmon"]
has_salmon = "salmon" in [a.lower() for a in aligners]

# 2) build all featureCounts inputs from the non-salmon list
aligner_fc_inputs = expand(
    cwd + "/results/{aligner}/all.{aligner}.{trimmer}_{pese}.fixcol2.featureCounts",
    aligner=non_salmon,
    trimmer=trimmers,
    pese=pese
)

# 3) optionally build the multi and multiFrac featureCounts
fctmulti_inputs = []
if runmulti.lower() == "true" and non_salmon:
    fctmulti_inputs = expand(
        cwd + "/results/{aligner}/all.{aligner}.{trimmer}_{pese}_multi.fixcol2.featureCounts",
        aligner=non_salmon,
        trimmer=trimmers,
        pese=pese
    )
if runmultifrac.lower() == "true" and non_salmon:
    fctmulti_inputs += expand(
        cwd + "/results/{aligner}/all.{aligner}.{trimmer}_{pese}_multifrac.fixcol2.featureCounts",
        aligner=non_salmon,
        trimmer=trimmers,
        pese=pese
    )

# 4) build salmon quant inputs if salmon is requested
salmon_inputs = []
if has_salmon:
    salmon_inputs = expand(
        "salmon/{trimmer}_{pese}/{unit.sample}.{unit.unit}/quant.sf",
        trimmer=trimmers,
        pese=pese,
        unit=units.itertuples()
    )

rule all:
    input:
        # pre-trim MultiQC
        expand("qc/multiqc_report_pretrim_{pese}.html", pese=pese),

        # post-trim MultiQC
        expand("qc/multiqc_report_{aligner}_{trimmer}_{pese}.html",
               aligner=aligners, trimmer=trimmers, pese=pese),

        # featureCounts for all non-salmon aligners
        aligner_fc_inputs,

        # salmon quant (if any)
        salmon_inputs,

        # multi & multifrac featureCounts
        fctmulti_inputs


include: "rules/qc.smk"
include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/featurecounts.smk"
