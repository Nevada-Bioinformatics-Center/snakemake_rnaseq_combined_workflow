rule featurecounts_onefile_star:
    input:
        sam=expand("star/{trimmer}/{unit.sample}.{unit.unit}/{unit.sample}.{unit.unit}_Aligned.sortedByCoord.out.bam", unit=units.itertuples(), trimmer=trimmers),
        annotation=config["ref"]["annotation"],
        fasta=config["ref"]["genomefa"]     # implicitly sets the -G flag
    output:
        multiext("results/star/all.star.{trimmer}",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts")
    log:
        "logs/star/{trimmer}/featurecount_all.log"
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra="{}".format(config["params"]["featurecounts"])
    threads: 16
    resources: time_min=220, mem_mb=20000, cpus=16
    wrapper:
        "0.73.0/bio/subread/featurecounts"

rule featurecounts_onefile_hisat2:
    input:
        sam=expand("hisat2/{trimmer}/{unit.sample}.{unit.unit}.sorted.bam", unit=units.itertuples(), trimmer=trimmers),
        annotation=config["ref"]["annotation"],
        fasta=config["ref"]["genomefa"]     # implicitly sets the -G flag
    output:
        multiext("results/hisat2/all.hisat2.{trimmer}",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts")
    log:
        "logs/hisat2/{trimmer}/featurecount_all.log"
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        #extra="-p -C -t exon -g gene_id --extraAttributes gene_name,transcript_id,transcript_name"   #-p to count fragments instead of reads, -C to not count fragments if mapped to different chrs, -t exon (only consider exons, -g look at the gene_id field to get the feature name
        extra="{}".format(config["params"]["featurecounts"])
    threads: 8
    resources: time_min=220, mem_mb=20000, cpus=8
    wrapper:
        "0.73.0/bio/subread/featurecounts"

rule fix_featurecounts:
    input:
	"results/{aligner}/featureCounts/all.featureCounts"
    output:
	"results/{aligner}/featureCounts/all.fixed.featureCounts"
    log:
        "logs/{aligner}/featurecount/fix.log"
    threads: 1
    resources: time_min=220, mem_mb=2000, cpus=1
    shell: 
        """
        python3 ../scripts/fix_featurecounts_output.py -f {input} -c 2 > {output} 2> {log}
        """


rule fix_featurecounts_hisat2:
    input:
        "results/hisat2/featureCounts/all.featureCounts"
    output:
	"results/hisat2/featureCounts/all.fixed.featureCounts"
    log:
        "logs/hisat2/featureCounts/fix.log"
    threads: 1
    resources: time_min=220, mem_mb=2000, cpus=1
    shell: 
        """
        python3 ../scripts/fix_featurecounts_output.py -f {input} -c 2 > {output} 2> {log}
        """

rule fix_featurecounts_star:
    input:
        "results/star/featureCounts/all.featureCounts"
    output:
	"results/star/featureCounts/all.fixed.featureCounts"
    log:
        "logs/star/featureCounts/fix.log"
    threads: 1
    resources: time_min=220, mem_mb=2000, cpus=1
    shell: 
        """
        python3 ../scripts/fix_featurecounts_output.py -f {input} -c 2 > {output} 2> {log}
        """
