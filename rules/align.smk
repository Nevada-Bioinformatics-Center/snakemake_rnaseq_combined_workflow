
def get_trim_fastq1(wildcards):
    fq1 = expand("trimmed/{trimmer}/{sample}.{unit}.1.fastq.gz", **wildcards)
    return fq1

def get_trim_fastq2(wildcards):
    fq2 = expand("trimmed/{trimmer}/{sample}.{unit}.2.fastq.gz", **wildcards)
    return fq2

##Star align
rule star_index:
    input:
        fasta = config["ref"]["genomefa"],
        gtf = config["ref"]["annotation"]
    output:
        directory(config["ref"]["index"] + "_star")
    params:
        #needed to change params below to build xenopus index
        #extra = "--limitGenomeGenerateRAM 60550893493 --genomeSAsparseD 3 --genomeSAindexNbases 12 -- genomeChrBinNbits 14"
        extra = ""
    threads: 16
    resources: time_min=480, mem_mb=40000, cpus=16
    log:
        "logs/star_index_genome.log"
    wrapper:
        #"0.71.1/bio/star/index"
        f"{wrappers_version}/bio/star/index"
            

rule star_align:
    input:
        fq1=get_trim_fastq1,
        fq2=get_trim_fastq2,
        genomedir=directory(config["ref"]["index"] + "_star")
    output:
        # see STAR manual for additional output files
        "star/{trimmer}/{sample}.{unit}/Aligned.sortedByCoord.out.bam",
        "star/{trimmer}/{sample}.{unit}/ReadsPerGene.out.tab"
    log:
        "logs/star/{trimmer}/{sample}.{unit}.log"
    params:
        # path to STAR reference genome index
        index=config["ref"]["index"] + "_star",
        # optional parameters
        extra="--outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --quantMode GeneCounts --sjdbGTFfile {} {}".format(
              config["ref"]["annotation"], config["params"]["star"])
    threads: 16
    resources: time_min=480, mem_mb=200000, cpus=16
    wrapper:
        #"0.73.0/bio/star/align"
        f"{wrappers_version}/bio/star/align"

rule symlink_bam:
    input:
        "star/{trimmer}/{sample}.{unit}/Aligned.sortedByCoord.out.bam" 
    output:
        "star/{trimmer}/{sample}.{unit}/{sample}.{unit}_Aligned.sortedByCoord.out.bam" 
    threads: 1
    shell:
        """
        ln -s Aligned.sortedByCoord.out.bam {output}
        """

rule samtools_index_star:
    input:
        "star/{trimmer}/{sample}.{unit}/{sample}.{unit}_Aligned.sortedByCoord.out.bam" 
    output:
        "star/{trimmer}/{sample}.{unit}/{sample}.{unit}_Aligned.sortedByCoord.out.bam.bai" 
    params:
        "" # optional params string
    resources: time_min=320, mem_mb=2000, cpus=1
    wrapper:
        #"0.73.0/bio/samtools/index"
        f"{wrappers_version}/bio/samtools/index"


##################
##Hisat2 align
##################
rule hisat2_extractexons:
    input:
        gtf = config["ref"]["annotation"]
    output:
        "hisat2_prep/genome.exons"
    log:
        "logs/hisat2_extract_exons.log"
    threads: 1
    resources: time_min=480, mem_mb=2000, cpus=1
    conda:
        "../envs/hisat2.yaml"
    shell:
        "hisat2_extract_exons.py {input.gtf} > {output} 2> {log}"

rule hisat2_extractsplicesites:
    input:
        gtf = config["ref"]["annotation"]
    output:
        "hisat2_prep/genome.ss"
    log:
        "logs/hisat2_extract_ss.log"
    threads: 1
    resources: time_min=480, mem_mb=2000, cpus=1
    conda:
        "../envs/hisat2.yaml"
    shell:
        "hisat2_extract_splice_sites.py {input.gtf} > {output} 2> {log}"

rule hisat2_index:
    input:
        fasta = config["ref"]["genomefa"],
        exons = "hisat2_prep/genome.exons",
        ss = "hisat2_prep/genome.ss"
    output:
        directory(config["ref"]["index"] + "_hisat2")
    params:
        prefix = config["ref"]["index"] + "_hisat2/genome"
    log:
        "logs/hisat2_index_genome.log"
    threads: 16
    resources: time_min=480, mem_mb=200000, cpus=16
    conda:
        "../envs/hisat2.yaml"
    shell:
        "mkdir {output} && hisat2-build -p {threads} --ss {input.ss} --exon {input.exons} {input.fasta} {params.prefix} 2> {log}"


rule hisat2_align:
    input:
        r1=get_trim_fastq1,
        r2=get_trim_fastq2,
        idx=config["ref"]["index"] + "_hisat2/"
    output:
        "hisat2/{trimmer}/{sample}.{unit}.bam"
    log:
        "logs/hisat2/{trimmer}/{sample}.{unit}.log"
    params:
      ## --new summary to allow multiqc parsing and 
      ## --dta to use XS BAM alignment information for stringtie downstream
        #extra="--new-summary --dta",
        extra="{}".format(config["params"]["hisat2"]),
        idx=config["ref"]["index"] + "_hisat2/genome",
    threads: 16
    wildcard_constraints:
        unit="rep\d+"
    resources: time_min=480, mem_mb=40000, cpus=16
    conda:
        "../envs/hisat2.yaml"
    shell:
        "(hisat2 --threads {threads} -x {params.idx} {params.extra} -1 {input.r1} -2 {input.r2} | samtools view -Sbh -o {output}) 2> {log}"

rule sambamba_sort:
    input:
        "hisat2/{trimmer}/{sample}.{unit}.bam"
    output:
        "hisat2/{trimmer}/{sample}.{unit}.sorted.bam"
    log:
        "logs/hisat2/{trimmer}/sambamba-sort/{sample}.{unit}.log"
    params: ""
    threads: 16 
    wildcard_constraints:
        unit="rep\d+"
    resources: time_min=480, mem_mb=20000, cpus=16
    wrapper:
        #"0.74.0/bio/sambamba/sort"
        f"{wrappers_version}/bio/sambamba/sort"

rule samtools_index_hisat2:
    input:
        "hisat2/{trimmer}/{sample}.{unit}.sorted.bam"
    output:
        "hisat2/{trimmer}/{sample}.{unit}.sorted.bam.bai"
    params:
        "" # optional params string
    resources: time_min=480, mem_mb=2000, cpus=1
    wrapper:
        #"0.73.0/bio/samtools/index"
        f"{wrappers_version}/bio/samtools/index"

