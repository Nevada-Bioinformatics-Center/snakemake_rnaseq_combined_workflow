##Optional rRNA removal (params: runsortmerna). Keeps non-rRNA reads for alignment.
def sortmerna_ref_args():
    return " ".join(f"--ref {r}" for r in sortmerna_refs)

rule sortmerna_pe:
    input:
        r1="trimmed/{trimmer}_pe/{sample}.{unit}.1.fastq.gz",
        r2="trimmed/{trimmer}_pe/{sample}.{unit}.2.fastq.gz",
        refs=sortmerna_refs,
    output:
        r1="sortmerna/{trimmer}_pe/{sample}.{unit}.1.fastq.gz",
        r2="sortmerna/{trimmer}_pe/{sample}.{unit}.2.fastq.gz",
        stats="qc/sortmerna/{trimmer}_pe/{sample}.{unit}.sortmerna.log",
    log:
        "logs/sortmerna/{trimmer}_pe/{sample}.{unit}.log"
    params:
        refs=sortmerna_ref_args(),
        extra=config["params"]["sortmerna"],
    threads: config["params"]["sortmernacpu"]
    resources: time_min=480, mem_mb=config["params"]["sortmernaram"], cpus=config["params"]["sortmernacpu"]
    conda:
        "../envs/sortmerna.yaml"
    shell:
        """
        wd=$(mktemp -d)
        trap 'rm -rf "$wd"' EXIT
        sortmerna {params.refs} --reads {input.r1} --reads {input.r2} \
            --workdir "$wd" --threads {threads} --fastx --out2 \
            --aligned "$wd/aligned" --other "$wd/other" {params.extra} > {log} 2>&1
        mv "$wd"/other_fwd.f* {output.r1}
        mv "$wd"/other_rev.f* {output.r2}
        mv "$wd/aligned.log" {output.stats}
        """

rule sortmerna_se:
    input:
        r1="trimmed/{trimmer}_se/{sample}.{unit}.1.fastq.gz",
        refs=sortmerna_refs,
    output:
        r1="sortmerna/{trimmer}_se/{sample}.{unit}.1.fastq.gz",
        stats="qc/sortmerna/{trimmer}_se/{sample}.{unit}.sortmerna.log",
    log:
        "logs/sortmerna/{trimmer}_se/{sample}.{unit}.log"
    params:
        refs=sortmerna_ref_args(),
        # --paired_in only applies to PE data
        extra=config["params"]["sortmerna"].replace("--paired_in", "").replace("--paired_out", ""),
    threads: config["params"]["sortmernacpu"]
    resources: time_min=480, mem_mb=config["params"]["sortmernaram"], cpus=config["params"]["sortmernacpu"]
    conda:
        "../envs/sortmerna.yaml"
    shell:
        """
        wd=$(mktemp -d)
        trap 'rm -rf "$wd"' EXIT
        sortmerna {params.refs} --reads {input.r1} \
            --workdir "$wd" --threads {threads} --fastx \
            --aligned "$wd/aligned" --other "$wd/other" {params.extra} > {log} 2>&1
        mv "$wd"/other.f* {output.r1}
        mv "$wd/aligned.log" {output.stats}
        """
