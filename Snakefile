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
    trimmer=r"[a-z]+",
    fct_mode = r"|_multi|_multifrac",

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
# A unit is paired-end when fq2 is filled in, otherwise single-end.
has_fq2 = units["fq2"].notna() & (units["fq2"] != "")
units["layout"] = has_fq2.map({True: "pe", False: "se"})
layouts_present = [l for l in ("pe", "se") if (units["layout"] == l).any()]

def units_for(layout):
    return units[units["layout"] == layout].itertuples()


aligners=config["params"]["aligners"].split(",")
trimmers=config["params"]["trimmers"].split(",")
runmulti=config["params"]["runfctmulti"]
runmultifrac=config["params"]["runfctmultifrac"]
runsortmerna = str(config["params"].get("runsortmerna", "false")).lower() == "true"
sortmerna_refs = [r.strip() for r in str(config["params"].get("sortmerna_refs", "")).split(",") if r.strip()]
if runsortmerna:
    if not sortmerna_refs:
        raise ValueError("runsortmerna is true but params: sortmerna_refs is empty in config.yaml")
    missing = [r for r in sortmerna_refs if not os.path.exists(r)]
    if missing:
        raise ValueError(f"sortmerna_refs not found: {missing}")
print("SortMeRNA:", runsortmerna)
runmosdepth = str(config["params"].get("runmosdepth", "false")).lower() == "true"
runrseqcgenebody = str(config["params"].get("runrseqcgenebody", "false")).lower() == "true"

def gtf_has_feature(path, feature):
    with open(path) as fh:
        for line in fh:
            if line[0] != "#" and line.split("\t", 3)[2:3] == [feature]:
                return True
    return False

# Bacterial annotations often have gene/CDS lines but no exon lines.
has_exons = gtf_has_feature(config["ref"]["annotation"], "exon")
feature_type = str(config["params"].get("feature_type", "auto"))
if feature_type.lower() == "auto":
    feature_type = "exon" if has_exons else "CDS"
if not gtf_has_feature(config["ref"]["annotation"], feature_type):
    raise ValueError(f"No '{feature_type}' features found in {config['ref']['annotation']}")

runpicard = str(config["params"].get("runpicardrnaseqmetrics", "auto")).lower()
if runpicard == "auto":
    runpicard = has_exons
else:
    runpicard = runpicard == "true"
    if runpicard and not has_exons:
        raise ValueError("runpicardrnaseqmetrics is true but the annotation has no exon features")
print("Coverage feature type:", feature_type, "| Picard RNA-seq metrics:", runpicard)
print("Aligners:", aligners)
print("Trimmers:", trimmers)
print("Layouts:", units["layout"].value_counts().to_dict())
cwd = os.getcwd() + "/"
print("Cwd:", cwd)

def strip_suffix(pattern, suffix):
    return pattern[: -len(suffix)]


wrappers_version="v7.0.0"


# build a list of which featureCounts variants to run
fct_modes = [""]  #
if runmulti.lower() == "true":
    fct_modes.append("_multi")
if runmultifrac.lower() == "true":
    fct_modes.append("_multifrac")



##Building rules to run
# 1) partition aligners
non_salmon = [a for a in aligners if a.lower() != "salmon"]
has_salmon = "salmon" in [a.lower() for a in aligners]


# 2) build salmon quant inputs if salmon is requested
salmon_inputs = []
if has_salmon:
    salmon_inputs = expand(
        "salmon/{trimmer}_{unit.layout}/{unit.sample}.{unit.unit}/quant.sf",
        trimmer=trimmers,
        unit=units.itertuples()
    )

rule all:
    input:
        # pre-trim MultiQC
        "qc/multiqc_report_pretrim.html",

        # post-trim MultiQC
        expand("qc/multiqc_report_{aligner}_{trimmer}.html",
               aligner=non_salmon, trimmer=trimmers),

        #  merged PE+SE featureCounts for all non salmon aligners (STAR + HISAT2),
        #    over every fct_mode ("", "_multi", "_multifrac")
        expand(
            cwd + "results/{aligner}/all.{aligner}.{trimmer}{fct_mode}.fixcol2.featureCounts",
            aligner=non_salmon,
            trimmer=trimmers,
            fct_mode=fct_modes
        ),

        # salmon quant (if any)
        salmon_inputs,


include: "rules/qc.smk"
include: "rules/trim.smk"
include: "rules/sortmerna.smk"
include: "rules/align.smk"
include: "rules/featurecounts.smk"
include: "rules/coverage.smk"
