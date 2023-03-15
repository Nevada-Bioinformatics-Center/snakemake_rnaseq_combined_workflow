rule featurecounts_onefile_star:
    input:
        samples=expand("star/{trimmer}/{unit.sample}.{unit.unit}/{unit.sample}.{unit.unit}_Aligned.sortedByCoord.out.bam", unit=units.itertuples(), trimmer=trimmers),
        annotation=config["ref"]["annotation"],
        fasta=config["ref"]["genomefa"]     # implicitly sets the -G flag
    output:
        multiext("results/star/all.star.{trimmer}",
                 ".featureCounts",
                 ".featureCounts.summary")
#                 ".featureCounts.jcounts")
    log:
        "logs/star/{trimmer}/featurecount_all.log"
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra="{}".format(config["params"]["featurecounts"])
    threads: 16
    resources: time_min=220, mem_mb=20000, cpus=16
    wrapper:
        f"{wrappers_version}/bio/subread/featurecounts"
        #"v1.2.0/bio/subread/featurecounts"
        #f"{wrappers_version}/bio/subread/featurecounts"

rule featurecounts_onefile_hisat2:
    input:
        samples=expand("hisat2/{trimmer}/{unit.sample}.{unit.unit}.sorted.bam", unit=units.itertuples(), trimmer=trimmers),
        annotation=config["ref"]["annotation"],
        fasta=config["ref"]["genomefa"]     # implicitly sets the -G flag
    output:
        multiext("results/hisat2/all.hisat2.{trimmer}",
                 ".featureCounts",
                 ".featureCounts.summary")
    log:
        "logs/hisat2/{trimmer}/featurecount_all.log"
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra="-P {}".format(config["params"]["featurecounts"])
    threads: 16
    resources: time_min=220, mem_mb=20000, cpus=16
    wrapper:
        f"{wrappers_version}/bio/subread/featurecounts"

rule featurecounts_onefile_hisat2_se:
    input:
        samples=expand("hisat2/{trimmer}_se/{unit.sample}.{unit.unit}.sorted.bam", unit=units.itertuples(), trimmer=trimmers),
        annotation=config["ref"]["annotation"],
        fasta=config["ref"]["genomefa"]     # implicitly sets the -G flag
    output:
        multiext("results/hisat2/all.hisat2.{trimmer}_se",
                 ".featureCounts",
                 ".featureCounts.summary")
    log:
        "logs/hisat2/{trimmer}_se/featurecount_all.log"
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra="{}".format(config["params"]["featurecountsse"])
    threads: 16
    resources: time_min=220, mem_mb=20000, cpus=16
    wrapper:
        f"{wrappers_version}/bio/subread/featurecounts"


rule fix_featurecounts_general:
    input:
        cwd+"/results/{aligner}/all.{aligner}.{trimmer}.featureCounts",
    output:
        cwd+"/results/{aligner}/all.{aligner}.{trimmer}.fixcol2.featureCounts",
    log:
        "logs/star/fct_fix_{aligner}_{trimmer}.log"
    threads: 1
    resources: time_min=220, mem_mb=2000, cpus=1
    shell: 
        "python3 scripts/fix_featurecounts_output.py -f {input} -c 2 > {output} 2> {log}"

rule fix_featurecounts_general_se:
    input:
        cwd+"/results/{aligner}/all.{aligner}.{trimmer}_se.featureCounts",
    output:
        cwd+"/results/{aligner}/all.{aligner}.{trimmer}_se.fixcol2.featureCounts",
    log:
        "logs/star/fct_fix_{aligner}_{trimmer}_se.log"
    threads: 1
    resources: time_min=220, mem_mb=2000, cpus=1
    shell: 
        "python3 scripts/fix_featurecounts_output.py -f {input} -c 2 > {output} 2> {log}"

rule fix_featurecounts_general_pe:
    input:
        cwd+"/results/{aligner}/all.{aligner}.{trimmer}_pe.featureCounts",
    output:
        cwd+"/results/{aligner}/all.{aligner}.{trimmer}_pe.fixcol2.featureCounts",
    log:
        "logs/star/fct_fix_{aligner}_{trimmer}_pe.log"
    threads: 1
    resources: time_min=220, mem_mb=2000, cpus=1
    shell: 
        "python3 scripts/fix_featurecounts_output.py -f {input} -c 2 > {output} 2> {log}"

#rule fix_featurecounts_star_fastp:
#    input:
#        "../results/star/all.star.fastp.featureCounts"
#    output:
#	"results/star/all.star.fastp.fixcol2.featureCounts"
#    log:
#        "logs/star/fct_fix_fastp.log"
#    threads: 1
#    resources: time_min=220, mem_mb=2000, cpus=1
#    shell: 
#        "python3 scripts/fix_featurecounts_output.py -f {input} -c 2 > {output} 2> {log}"
#
#rule fix_featurecounts_hisat2_fastp:
#    input:
#        "../results/hisat2/all.hisat2.fastp.featureCounts"
#    output:
#	"results/hisat2/all.hisat2.fastp.fixcol2.featureCounts"
#    log:
#        "logs/hisat2/fct_fix_fastp.log"
#    threads: 1
#    resources: time_min=220, mem_mb=2000, cpus=1
#    shell: 
#        "python3 ../scripts/fix_featurecounts_output.py -f {input} -c 2 > {output} 2> {log}"
