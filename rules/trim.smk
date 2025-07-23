def get_fastq(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

def get_fastq1(wildcards):
    if pese == "se":
        return [(units.loc[(wildcards.sample, wildcards.unit), ["fq1"]].dropna().item())]
    else:
        return (units.loc[(wildcards.sample, wildcards.unit), ["fq1"]].dropna().item())

def get_fastq2(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), ["fq2"]].dropna().item()

rule trimmomatic_pe:
    input:
        r1=get_fastq1,
        r2=get_fastq2
    output:
        r1="trimmed/trimmomatic_pe/{sample}.{unit}.1.fastq.gz",
        r2="trimmed/trimmomatic_pe/{sample}.{unit}.2.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired="trimmed/trimmomatic_pe/{sample}.{unit}.1.unpaired.fastq.gz",
        r2_unpaired="trimmed/trimmomatic_pe/{sample}.{unit}.2.unpaired.fastq.gz"
    log:
        "logs/trimmomatic_pe/{sample}.{unit}.log"
    params:
        # list of trimmers (see manual)
        trimmer = [f"ILLUMINACLIP:{config['ref']['adapter']}:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36"],
        extra="",
        compression_level="-9"
    threads: config["params"]["trimmomaticcpu"]
    resources: time_min=480, mem_mb=config["params"]["trimmomaticram"], cpus=config["params"]["trimmomaticcpu"]
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    wrapper:
        f"{wrappers_version}/bio/trimmomatic/pe"

rule trimmomatic_se:
    input:
        r1=get_fastq1
    output:
        r1="trimmed/trimmomatic_se/{sample}.{unit}.1.fastq.gz",
        r1_unpaired="trimmed/trimmomatic_se/{sample}.{unit}.1.unpaired.fastq.gz",
    log:
        "logs/trimmomatic_se/{sample}.{unit}.log"
    params:
        # list of trimmers (see manual)
        trimmer = [f"ILLUMINACLIP:{config['ref']['adapter']}:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36"],
        extra="",
        compression_level="-9"
    threads: config["params"]["trimmomaticcpu"]
    resources: time_min=480, mem_mb=config["params"]["trimmomaticram"], cpus=config["params"]["trimmomaticcpu"]
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    wrapper:
        f"{wrappers_version}/bio/trimmomatic/se"

rule fastp_pe:
    input:
        sample=get_fastq
    output:
        trimmed=["trimmed/fastp_pe/{sample}.{unit}.1.fastq.gz", "trimmed/fastp_pe/{sample}.{unit}.2.fastq.gz"],
        html="report/fastp_pe/{sample}.{unit}.html",
        json="report/fastp_pe/{sample}.{unit}.fastp.json"
    log:
        "logs/fastp_pe/{sample}.{unit}.log"
    params:
        #adapters="--adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        adapters="--detect_adapter_for_pe",
        extra=config["params"]["fastp"]
    threads: config["params"]["fastpcpu"]
    resources: time_min=480, mem_mb=config["params"]["fastpram"], cpus=config["params"]["fastpcpu"]
    wrapper:
        f"{wrappers_version}/bio/fastp"

rule fastp_se:
    input:
        sample=get_fastq1
    output:
        trimmed="trimmed/fastp_se/{sample}.{unit}.1.fastq.gz",
        failed="trimmed/fastp_se/{sample}.{unit}.failed.fastq.gz",
        html="report/fastp_se/{sample}.{unit}.html",
        json="report/fastp_se/{sample}.{unit}.fastp.json"
    log:
        "logs/fastp_se/{sample}.{unit}.log"
    params:
        adapters="",
        extra=config["params"]["fastp"]
    threads: config["params"]["fastpcpu"]
    resources: time_min=480, mem_mb=config["params"]["fastpram"], cpus=config["params"]["fastpcpu"]
    wrapper:
        f"{wrappers_version}/bio/fastp"

rule trim_galore_pe:
    input:
        get_fastq
    output:
        fasta_fwd="trimmed/trimgalore_pe/{sample}.{unit}.1.fastq.gz",
        report_fwd="trimmed/trimgalore_pe/{sample}.{unit}.1_trimming_report.txt",
        fasta_rev="trimmed/trimgalore_pe/{sample}.{unit}.2.fastq.gz",
        report_rev="trimmed/trimgalore_pe/{sample}.{unit}.2_trimming_report.txt",
    params:
        extra="--illumina -q 20",
    log:
        "logs/trimgalore/{sample}.{unit}.log",
    threads: config["params"]["trimgalorecpu"]
    resources: time_min=480, mem_mb=config["params"]["trimgaloreram"], cpus=config["params"]["trimgalorecpu"]
    wrapper:
        f"{wrappers_version}/bio/trim_galore/pe"

rule trim_galore_se:
    input:
        get_fastq
    output:
        fasta="trimmed/trimgalore_se/{sample}.{unit}.1.fastq.gz",
        report="trimmed/trimgalore_se/{sample}.{unit}.1_trimming_report.txt",
    params:
        extra="--illumina -q 20",
    log:
        "logs/trimgalore/{sample}.{unit}.log",
    threads: config["params"]["trimgalorecpu"]
    resources: time_min=480, mem_mb=config["params"]["trimgaloreram"], cpus=config["params"]["trimgalorecpu"]
    wrapper:
        f"{wrappers_version}/bio/trim_galore/se"
