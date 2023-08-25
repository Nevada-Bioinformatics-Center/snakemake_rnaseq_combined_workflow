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

rule featurecounts_onefile_star_se:
    input:
        samples=expand("star/{trimmer}_se/{unit.sample}.{unit.unit}/{unit.sample}.{unit.unit}_Aligned.sortedByCoord.out.bam", unit=units.itertuples(), trimmer=trimmers),
        annotation=config["ref"]["annotation"],
        fasta=config["ref"]["genomefa"]     # implicitly sets the -G flag
    output:
        multiext("results/star/all.star.{trimmer}_se",
                 ".featureCounts",
                 ".featureCounts.summary")
#                 ".featureCounts.jcounts")
    log:
        "logs/star/{trimmer}/featurecount_all_se.log"
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra="{}".format(config["params"]["featurecountsse"])
    threads: 16
    resources: time_min=220, mem_mb=20000, cpus=16
    wrapper:
        f"{wrappers_version}/bio/subread/featurecounts"

rule featurecounts_onefile_star_pe:
    input:
        samples=expand("star/{trimmer}_pe/{unit.sample}.{unit.unit}/{unit.sample}.{unit.unit}_Aligned.sortedByCoord.out.bam", unit=units.itertuples(), trimmer=trimmers),
        annotation=config["ref"]["annotation"],
        fasta=config["ref"]["genomefa"]     # implicitly sets the -G flag
    output:
        multiext("results/star/all.star.{trimmer}_pe",
                 ".featureCounts",
                 ".featureCounts.summary")
#                 ".featureCounts.jcounts")
    log:
        "logs/star/{trimmer}/featurecount_all_pe.log"
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra="{}".format(config["params"]["featurecounts"])
    threads: 16
    resources: time_min=220, mem_mb=20000, cpus=16
    wrapper:
        f"{wrappers_version}/bio/subread/featurecounts"

rule featurecounts_onefile_star_pe_multi:
    input:
        samples=expand("star/{trimmer}_pe/{unit.sample}.{unit.unit}/{unit.sample}.{unit.unit}_Aligned.sortedByCoord.out.bam", unit=units.itertuples(), trimmer=trimmers),
        annotation=config["ref"]["annotation"],
        fasta=config["ref"]["genomefa"]     # implicitly sets the -G flag
    output:
        multiext("results/star/all.star.{trimmer}_pe_multi",
                 ".featureCounts",
                 ".featureCounts.summary")
#                 ".featureCounts.jcounts")
    log:
        "logs/star/{trimmer}/featurecount_all_pe_multi.log"
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra="{}".format(config["params"]["featurecountsmulti"])
    threads: 16
    resources: time_min=220, mem_mb=20000, cpus=16
    wrapper:
        f"{wrappers_version}/bio/subread/featurecounts"

rule featurecounts_onefile_star_pe_multifrac:
    input:
        samples=expand("star/{trimmer}_pe/{unit.sample}.{unit.unit}/{unit.sample}.{unit.unit}_Aligned.sortedByCoord.out.bam", unit=units.itertuples(), trimmer=trimmers),
        annotation=config["ref"]["annotation"],
        fasta=config["ref"]["genomefa"]     # implicitly sets the -G flag
    output:
        multiext("results/star/all.star.{trimmer}_pe_multifrac",
                 ".featureCounts",
                 ".featureCounts.summary")
#                 ".featureCounts.jcounts")
    log:
        "logs/star/{trimmer}/featurecount_all_pe_multifrac.log"
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra="{}".format(config["params"]["featurecountsmultifrac"])
    threads: 16
    resources: time_min=220, mem_mb=20000, cpus=16
    wrapper:
        f"{wrappers_version}/bio/subread/featurecounts"

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

rule featurecounts_onefile_hisat2_pe:
    input:
        samples=expand("hisat2/{trimmer}_pe/{unit.sample}.{unit.unit}.sorted.bam", unit=units.itertuples(), trimmer=trimmers),
        annotation=config["ref"]["annotation"],
        fasta=config["ref"]["genomefa"]     # implicitly sets the -G flag
    output:
        multiext("results/hisat2/all.hisat2.{trimmer}_pe",
                 ".featureCounts",
                 ".featureCounts.summary")
    log:
        "logs/hisat2/{trimmer}_pe/featurecount_all.log"
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra="{}".format(config["params"]["featurecounts"])
    threads: 16
    resources: time_min=220, mem_mb=20000, cpus=16
    wrapper:
        f"{wrappers_version}/bio/subread/featurecounts"

rule featurecounts_onefile_hisat2_pe_multi:
    input:
        samples=expand("hisat2/{trimmer}_pe/{unit.sample}.{unit.unit}.sorted.bam", unit=units.itertuples(), trimmer=trimmers),
        annotation=config["ref"]["annotation"],
        fasta=config["ref"]["genomefa"]     # implicitly sets the -G flag
    output:
        multiext("results/hisat2/all.hisat2.{trimmer}_pe_multi",
                 ".featureCounts",
                 ".featureCounts.summary")
    log:
        "logs/hisat2/{trimmer}_pe/featurecount_all_multi.log"
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra="{}".format(config["params"]["featurecountsmulti"])
    threads: 16
    resources: time_min=220, mem_mb=20000, cpus=16
    wrapper:
        f"{wrappers_version}/bio/subread/featurecounts"

rule featurecounts_onefile_hisat2_pe_multifrac:
    input:
        samples=expand("hisat2/{trimmer}_pe/{unit.sample}.{unit.unit}.sorted.bam", unit=units.itertuples(), trimmer=trimmers),
        annotation=config["ref"]["annotation"],
        fasta=config["ref"]["genomefa"]     # implicitly sets the -G flag
    output:
        multiext("results/hisat2/all.hisat2.{trimmer}_pe_multifrac",
                 ".featureCounts",
                 ".featureCounts.summary")
    log:
        "logs/hisat2/{trimmer}_pe/featurecount_all_multifrac.log"
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra="{}".format(config["params"]["featurecountsmultifrac"])
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
        "results/{aligner}/all.{aligner}.{trimmer}_se.featureCounts",
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
        "results/{aligner}/all.{aligner}.{trimmer}_pe.featureCounts",
    output:
        cwd+"/results/{aligner}/all.{aligner}.{trimmer}_pe.fixcol2.featureCounts",
    log:
        "logs/star/fct_fix_{aligner}_{trimmer}_pe.log"
    threads: 1
    resources: time_min=220, mem_mb=2000, cpus=1
    shell: 
        "python3 scripts/fix_featurecounts_output.py -f {input} -c 2 > {output} 2> {log}"

rule fix_featurecounts_general_pe_multi:
    input:
        "results/{aligner}/all.{aligner}.{trimmer}_pe_multi.featureCounts",
    output:
        cwd+"/results/{aligner}/all.{aligner}.{trimmer}_pe_multi.fixcol2.featureCounts",
    log:
        "logs/star/fct_fix_{aligner}_{trimmer}_pe_multi.log"
    threads: 1
    resources: time_min=220, mem_mb=2000, cpus=1
    shell: 
        "python3 scripts/fix_featurecounts_output.py -f {input} -c 2 > {output} 2> {log}"

rule fix_featurecounts_general_pe_multifrac:
    input:
        "results/{aligner}/all.{aligner}.{trimmer}_pe_multifrac.featureCounts",
    output:
        cwd+"/results/{aligner}/all.{aligner}.{trimmer}_pe_multifrac.fixcol2.featureCounts",
    log:
        "logs/star/fct_fix_{aligner}_{trimmer}_pe_multifrac.log"
    threads: 1
    resources: time_min=220, mem_mb=2000, cpus=1
    shell: 
        "python3 scripts/fix_featurecounts_output.py -f {input} -c 2 > {output} 2> {log}"
