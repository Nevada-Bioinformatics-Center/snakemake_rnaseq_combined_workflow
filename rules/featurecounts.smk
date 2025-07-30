rule featurecounts_onefile_star:
    input:
        samples=lambda wildcards: expand("star/{trimmer}_{pese}/{unit.sample}.{unit.unit}/{unit.sample}.{unit.unit}_Aligned.sortedByCoord.out.bam", unit=units.itertuples(), trimmer=wildcards.trimmer, pese=wildcards.pese),
        annotation=config["ref"]["annotation"],
        fasta=config["ref"]["genomefa"]     # implicitly sets the -G flag
    output:
        multiext("results/star/all.star.{trimmer}_{pese}",
                 ".featureCounts",
                 ".featureCounts.summary")
    log:
        "logs/star/{trimmer}/featurecount_all_{pese}.log"
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra=lambda wc: (
            # for PE data, add -p -B -C
            "-p -B -C " + config["params"]["featurecounts"]
            if wc.pese == "pe"
            # for SE data, just use the base params
            else config["params"]["featurecounts"]
        ),
    threads: config["params"]["featurecountscpu"]
    resources: time_min=220, mem_mb=config["params"]["featurecountsram"], cpus=config["params"]["featurecountscpu"]
    wrapper:
        f"{wrappers_version}/bio/subread/featurecounts"

rule featurecounts_onefile_star_multi:
    input:
        samples=lambda wildcards: expand("star/{trimmer}_{pese}/{unit.sample}.{unit.unit}/{unit.sample}.{unit.unit}_Aligned.sortedByCoord.out.bam", unit=units.itertuples(), trimmer=wildcards.trimmer, pese=wildcards.pese),
        annotation=config["ref"]["annotation"],
        fasta=config["ref"]["genomefa"]     # implicitly sets the -G flag
    output:
        multiext("results/star/all.star.{trimmer}_{pese}_multi",
                 ".featureCounts",
                 ".featureCounts.summary")
    log:
        "logs/star/{trimmer}/featurecount_all_{pese}_multi.log"
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra=lambda wc: (
            # for PE data, add -p -B -C
            "-p -B -C " + config["params"]["featurecountsmulti"]
            if wc.pese == "pe"
            # for SE data, just use the base params
            else config["params"]["featurecountsmulti"]
        ),
    threads: config["params"]["featurecountscpu"]
    resources: time_min=220, mem_mb=config["params"]["featurecountsram"], cpus=config["params"]["featurecountscpu"]
    wrapper:
        f"{wrappers_version}/bio/subread/featurecounts"

rule featurecounts_onefile_star_multifrac:
    input:
        samples=lambda wildcards: expand("star/{trimmer}_{pese}/{unit.sample}.{unit.unit}/{unit.sample}.{unit.unit}_Aligned.sortedByCoord.out.bam", unit=units.itertuples(), trimmer=wildcards.trimmer, pese=wildcards.pese),
        annotation=config["ref"]["annotation"],
        fasta=config["ref"]["genomefa"]     # implicitly sets the -G flag
    output:
        multiext("results/star/all.star.{trimmer}_{pese}_multifrac",
                 ".featureCounts",
                 ".featureCounts.summary")
    log:
        "logs/star/{trimmer}/featurecount_all_{pese}_multifrac.log"
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra=lambda wc: (
            # for PE data, add -p -B -C
            "-p -B -C " + config["params"]["featurecountsmultifrac"]
            if wc.pese == "pe"
            # for SE data, just use the base params
            else config["params"]["featurecountsmultifrac"]
        ),
    threads: config["params"]["featurecountscpu"]
    resources: time_min=220, mem_mb=config["params"]["featurecountsram"], cpus=config["params"]["featurecountscpu"]
    wrapper:
        f"{wrappers_version}/bio/subread/featurecounts"



rule featurecounts_onefile_hisat2:
    input:
        samples=lambda wildcards: expand("hisat2/{trimmer}_{pese}/{unit.sample}.{unit.unit}.sorted.bam", unit=units.itertuples(), trimmer=wildcards.trimmer, pese=wildcards.pese),
        annotation=config["ref"]["annotation"],
        fasta=config["ref"]["genomefa"]
    output:
        multiext("results/hisat2/all.hisat2.{trimmer}_{pese}",
                 ".featureCounts",
                 ".featureCounts.summary")
    log:
        "logs/hisat2/{trimmer}_{pese}/featurecount_all.log"
    params:
        tmp_dir="",
        r_path="",
        extra=lambda wc: (
            # for PE data, add -p -B -C
            "-p -B -C " + config["params"]["featurecounts"]
            if wc.pese == "pe"
            # for SE data, just use the base params
            else config["params"]["featurecounts"]
        ),
    threads: config["params"]["featurecountscpu"]
    resources: time_min=220, mem_mb=config["params"]["featurecountsram"], cpus=config["params"]["featurecountscpu"]
    wrapper:
        f"{wrappers_version}/bio/subread/featurecounts"


rule featurecounts_onefile_hisat2_multi:
    input:
        samples=lambda wildcards: expand("hisat2/{trimmer}_{pese}/{unit.sample}.{unit.unit}.sorted.bam", unit=units.itertuples(), trimmer=wildcards.trimmer, pese=wildcards.pese),
        annotation=config["ref"]["annotation"],
        fasta=config["ref"]["genomefa"]     # implicitly sets the -G flag
    output:
        multiext("results/hisat2/all.hisat2.{trimmer}_{pese}_multi",
                 ".featureCounts",
                 ".featureCounts.summary")
    log:
        "logs/hisat2/{trimmer}_{pese}/featurecount_all_multi.log"
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra=lambda wc: (
            # for PE data, add -p -B -C
            "-p -B -C " + config["params"]["featurecountsmulti"]
            if wc.pese == "pe"
            # for SE data, just use the base params
            else config["params"]["featurecountsmulti"]
        ),
    threads: config["params"]["featurecountscpu"]
    resources: time_min=220, mem_mb=config["params"]["featurecountsram"], cpus=config["params"]["featurecountscpu"]
    wrapper:
        f"{wrappers_version}/bio/subread/featurecounts"

rule featurecounts_onefile_hisat2_multifrac:
    input:
        samples=lambda wildcards: expand("hisat2/{trimmer}_{pese}/{unit.sample}.{unit.unit}.sorted.bam", unit=units.itertuples(), trimmer=wildcards.trimmer, pese=wildcards.pese),
        annotation=config["ref"]["annotation"],
        fasta=config["ref"]["genomefa"]     # implicitly sets the -G flag
    output:
        multiext("results/hisat2/all.hisat2.{trimmer}_{pese}_multifrac",
                 ".featureCounts",
                 ".featureCounts.summary")
    log:
        "logs/hisat2/{trimmer}_{pese}/featurecount_all_multifrac.log"
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra=lambda wc: (
            # for PE data, add -p -B -C
            "-p -B -C " + config["params"]["featurecountsmultifrac"]
            if wc.pese == "pe"
            # for SE data, just use the base params
            else config["params"]["featurecountsmultifrac"]
        ),
    threads: config["params"]["featurecountscpu"]
    resources: time_min=220, mem_mb=config["params"]["featurecountsram"], cpus=config["params"]["featurecountscpu"]
    wrapper:
        f"{wrappers_version}/bio/subread/featurecounts"


rule fix_featurecounts_general:
    input:
        "results/{aligner}/all.{aligner}.{trimmer}_{pese}.featureCounts",
    output:
        cwd+"/results/{aligner}/all.{aligner}.{trimmer}_{pese}.fixcol2.featureCounts",
    log:
        "logs/star/fct_fix_{aligner}_{trimmer}_{pese}.log"
    threads: 1
    resources: time_min=220, mem_mb=8000, cpus=1
    shell: 
        "python3 scripts/fix_featurecounts_output.py -f {input} -c 2 > {output} 2> {log}"


rule fix_featurecounts_general_multi:
    input:
        "results/{aligner}/all.{aligner}.{trimmer}_{pese}_multi.featureCounts",
    output:
        cwd+"/results/{aligner}/all.{aligner}.{trimmer}_{pese}_multi.fixcol2.featureCounts",
    log:
        "logs/star/fct_fix_{aligner}_{trimmer}_{pese}_multi.log"
    threads: 1
    resources: time_min=220, mem_mb=8000, cpus=1
    shell: 
        "python3 scripts/fix_featurecounts_output.py -f {input} -c 2 > {output} 2> {log}"

rule fix_featurecounts_general_multifrac:
    input:
        "results/{aligner}/all.{aligner}.{trimmer}_{pese}_multifrac.featureCounts",
    output:
        cwd+"/results/{aligner}/all.{aligner}.{trimmer}_{pese}_multifrac.fixcol2.featureCounts",
    log:
        "logs/star/fct_fix_{aligner}_{trimmer}_{pese}_multifrac.log"
    threads: 1
    resources: time_min=220, mem_mb=8000, cpus=1
    shell: 
        "python3 scripts/fix_featurecounts_output.py -f {input} -c 2 > {output} 2> {log}"

