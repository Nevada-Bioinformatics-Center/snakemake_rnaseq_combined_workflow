##Coverage QC on aligned BAMs (params: runmosdepth, runrseqcgenebody)
rule exons_bed:
    input:
        gtf="results/ref/prepared_annotation.gtf"
    output:
        bed="results/ref/exons.merged.bed"
    log:
        "logs/refs/exons_bed.log"
    conda:
        "../envs/mosdepth.yaml"
    shell:
        """
        awk 'BEGIN {{OFS="\\t"}} $3=="exon" {{print $1, $4-1, $5}}' {input.gtf} \
            | sort -k1,1 -k2,2n | bedtools merge -i - > {output.bed} 2> {log}
        """

rule gtf_to_bed12:
    input:
        "results/ref/prepared_annotation.gtf"
    output:
        db="results/ref/prepared_annotation.gffutils.db",
        bed="results/ref/transcripts.bed12"
    log:
        "logs/refs/gtf_to_bed12.log"
    conda:
        "../envs/gffutils.yaml"
    script:
        "../scripts/gtf2bed.py"

rule mosdepth:
    input:
        bam=get_picard_input,
        bai=lambda wc: get_picard_input(wc) + ".bai",
        regions="results/ref/exons.merged.bed"
    output:
        summary="qc/mosdepth/{aligner}/{trimmer}_{pese}/{sample}.{unit}.mosdepth.summary.txt",
        global_dist="qc/mosdepth/{aligner}/{trimmer}_{pese}/{sample}.{unit}.mosdepth.global.dist.txt",
        region_dist="qc/mosdepth/{aligner}/{trimmer}_{pese}/{sample}.{unit}.mosdepth.region.dist.txt"
    wildcard_constraints:
        aligner=r"star|hisat2"
    log:
        "logs/mosdepth/{aligner}/{trimmer}_{pese}/{sample}.{unit}.log"
    params:
        prefix=lambda wc, output: output.summary[: -len(".mosdepth.summary.txt")],
        extra=config["params"]["mosdepth"]
    threads: config["params"]["mosdepthcpu"]
    resources: runtime=320, mem_mb=config["params"]["mosdepthram"], cpus=config["params"]["mosdepthcpu"]
    conda:
        "../envs/mosdepth.yaml"
    shell:
        "mosdepth --threads {threads} --no-per-base --by {input.regions} {params.extra} "
        "{params.prefix} {input.bam} > {log} 2>&1"

rule rseqc_genebody_coverage:
    input:
        bam=get_picard_input,
        bai=lambda wc: get_picard_input(wc) + ".bai",
        bed="results/ref/transcripts.bed12"
    output:
        "qc/rseqc/{aligner}/{trimmer}_{pese}/{sample}.{unit}.geneBodyCoverage.txt"
    wildcard_constraints:
        aligner=r"star|hisat2"
    log:
        "logs/rseqc/genebody/{aligner}/{trimmer}_{pese}/{sample}.{unit}.log"
    params:
        prefix=lambda wc, output: output[0][: -len(".geneBodyCoverage.txt")]
    threads: 1
    resources: runtime=480, mem_mb=config["params"]["rseqcram"], cpus=1
    conda:
        "../envs/rseqc.yaml"
    shell:
        "geneBody_coverage.py -i {input.bam} -r {input.bed} -o {params.prefix} > {log} 2>&1"
