
def get_trim_fastq1(wildcards):
    if pese == "pe":
        fq1 = expand("trimmed/{trimmer}_pe/{sample}.{unit}.1.fastq.gz", **wildcards)
    else:
        fq1 = expand("trimmed/{trimmer}_se/{sample}.{unit}.1.fastq.gz", **wildcards)
    return fq1

def get_trim_fastq2(wildcards):
    fq2 = expand("trimmed/{trimmer}_pe/{sample}.{unit}.2.fastq.gz", **wildcards)
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
        #extra = ""
        extra = config["params"]["starindex"],
    threads: 16
    resources: time_min=480, mem_mb=129000, cpus=16
    log:
        "logs/star_index_genome.log"
    wrapper:
        #"0.71.1/bio/star/index"
        f"{wrappers_version}/bio/star/index"
            

rule star_align_pe:
    input:
        fq1=get_trim_fastq1,
        fq2=get_trim_fastq2,
        genomedir=directory(config["ref"]["index"] + "_star")
    output:
        # see STAR manual for additional output files
        aln="star/{trimmer}_pe/{sample}.{unit}/Aligned.sortedByCoord.out.bam",
        log="star/{trimmer}_pe/{sample}.{unit}/Log.out",
        sj="star/{trimmer}_pe/{sample}.{unit}/ReadsPerGene.out.tab"
    log:
        "logs/star/{trimmer}_pe/{sample}.{unit}.log"
    params:
        # path to STAR reference genome index
        idx=config["ref"]["index"] + "_star",
        # optional parameters
        extra="--outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --quantMode GeneCounts --sjdbGTFfile {} {}".format(
              config["ref"]["annotation"], config["params"]["star"])
    threads: 16
    resources: time_min=480, mem_mb=50000, cpus=16
    wrapper:
        f"{wrappers_version}/bio/star/align"

rule symlink_bam_pe:
    input:
        "star/{trimmer}_pe/{sample}.{unit}/Aligned.sortedByCoord.out.bam" 
    output:
        "star/{trimmer}_pe/{sample}.{unit}/{sample}.{unit}_Aligned.sortedByCoord.out.bam" 
    threads: 1
    shell:
        """
        ln -s Aligned.sortedByCoord.out.bam {output}
        """

rule samtools_index_star_pe:
    input:
        "star/{trimmer}_pe/{sample}.{unit}/{sample}.{unit}_Aligned.sortedByCoord.out.bam" 
    output:
        "star/{trimmer}_pe/{sample}.{unit}/{sample}.{unit}_Aligned.sortedByCoord.out.bam.bai" 
    params:
        "" # optional params string
    resources: time_min=320, mem_mb=2000, cpus=1
    wrapper:
        f"{wrappers_version}/bio/samtools/index"

rule star_align_se:
    input:
        fq1=get_trim_fastq1,
        genomedir=directory(config["ref"]["index"] + "_star")
    output:
        # see STAR manual for additional output files
        aln="star/{trimmer}_se/{sample}.{unit}/Aligned.sortedByCoord.out.bam",
        log="star/{trimmer}_se/{sample}.{unit}/Log.out",
        sj="star/{trimmer}_se/{sample}.{unit}/ReadsPerGene.out.tab"
    log:
        "logs/star/{trimmer}/{sample}.{unit}.log"
    params:
        idx=config["ref"]["index"] + "_star",
        extra="--outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --quantMode GeneCounts --sjdbGTFfile {} {}".format(
              config["ref"]["annotation"], config["params"]["star"])
    threads: 16
    resources: time_min=480, mem_mb=200000, cpus=16
    wrapper:
        f"{wrappers_version}/bio/star/align"

rule symlink_bam_se:
    input:
        "star/{trimmer}_se/{sample}.{unit}/Aligned.sortedByCoord.out.bam" 
    output:
        "star/{trimmer}_se/{sample}.{unit}/{sample}.{unit}_Aligned.sortedByCoord.out.bam" 
    threads: 1
    shell:
        """
        ln -s Aligned.sortedByCoord.out.bam {output}
        """

rule samtools_index_star_se:
    input:
        "star/{trimmer}_se/{sample}.{unit}/{sample}.{unit}_Aligned.sortedByCoord.out.bam" 
    output:
        "star/{trimmer}_se/{sample}.{unit}/{sample}.{unit}_Aligned.sortedByCoord.out.bam.bai" 
    params:
        "" # optional params string
    resources: time_min=320, mem_mb=2000, cpus=1
    wrapper:
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

rule hisat2_index_noexons:
    input:
        fasta = config["ref"]["genomefa"],
    output:
        directory(config["ref"]["index"] + "_hisat2"),
        config["ref"]["index"] + "_hisat2/genome.1.ht2"
    params:
        prefix = config["ref"]["index"] + "_hisat2/genome",
        odir=config["ref"]["index"] + "_hisat2"
    log:
        "logs/hisat2_index_genome.log"
    threads: 16
    resources: time_min=480, mem_mb=200000, cpus=16
    conda:
        "../envs/hisat2.yaml"
    shell:
        #"mkdir {params.odir} && hisat2-build -p {threads} {input.fasta} {params.prefix} 2> {log}"
        "hisat2-build -p {threads} {input.fasta} {params.prefix} 2> {log}"

#rule hisat2_index:
#    input:
#        fasta = config["ref"]["genomefa"],
#        exons = "hisat2_prep/genome.exons",
#        ss = "hisat2_prep/genome.ss"
#    output:
#        directory(config["ref"]["index"] + "_hisat2")
#    params:
#        prefix = config["ref"]["index"] + "_hisat2/genome"
#    log:
#        "logs/hisat2_index_genome.log"
#    threads: 16
#    resources: time_min=480, mem_mb=200000, cpus=16
#    conda:
#        "../envs/hisat2.yaml"
#    shell:
#        "mkdir {output} && hisat2-build -p {threads} --ss {input.ss} --exon {input.exons} {input.fasta} {params.prefix} 2> {log}"


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
        unit=r"rep\d+"
    resources: time_min=480, mem_mb=40000, cpus=16
    conda:
        "../envs/hisat2.yaml"
    shell:
        "(hisat2 --threads {threads} -x {params.idx} {params.extra} -1 {input.r1} -2 {input.r2} | samtools view -Sbh -o {output}) 2> {log}"

rule hisat2_align_pe:
    input:
        r1=get_trim_fastq1,
        r2=get_trim_fastq2,
        idx=config["ref"]["index"] + "_hisat2/"
    output:
        "hisat2/{trimmer}_pe/{sample}.{unit}.bam"
    log:
        "logs/hisat2/{trimmer}_pe/{sample}.{unit}.log"
    params:
      ## --new summary to allow multiqc parsing and 
      ## --dta to use XS BAM alignment information for stringtie downstream
        #extra="--new-summary --dta",
        extra="{}".format(config["params"]["hisat2"]),
        idx=config["ref"]["index"] + "_hisat2/genome",
    threads: 16
    wildcard_constraints:
        unit=r"rep\d+"
    resources: time_min=480, mem_mb=40000, cpus=16
    conda:
        "../envs/hisat2.yaml"
    shell:
        "(hisat2 --threads {threads} -x {params.idx} {params.extra} -1 {input.r1} -2 {input.r2} | samtools view -Sbh -o {output}) 2> {log}"

rule hisat2_align_se:
    input:
        r1=get_trim_fastq1,
        idx=config["ref"]["index"] + "_hisat2/"
    output:
        "hisat2/{trimmer}_se/{sample}.{unit}.bam"
    log:
        "logs/hisat2/{trimmer}_se/{sample}.{unit}.log"
    params:
      ## --new summary to allow multiqc parsing and 
      ## --dta to use XS BAM alignment information for stringtie downstream
        #extra="--new-summary --dta",
        extra="{}".format(config["params"]["hisat2"]),
        idx=config["ref"]["index"] + "_hisat2/genome",
    threads: 16
    wildcard_constraints:
        unit=r"rep\d+"
    resources: time_min=480, mem_mb=40000, cpus=16
    conda:
        "../envs/hisat2.yaml"
    shell:
        "(hisat2 --threads {threads} -x {params.idx} {params.extra} -U {input.r1} | samtools view -Sbh -o {output}) 2> {log}"

rule sambamba_sort_se:
    input:
        "hisat2/{trimmer}_se/{sample}.{unit}.bam"
    output:
        "hisat2/{trimmer}_se/{sample}.{unit}.sorted.bam"
    log:
        "logs/hisat2/{trimmer}_se/sambamba-sort/{sample}.{unit}.log"
    params: ""
    threads: 16 
    wildcard_constraints:
        unit=r"rep\d+"
    resources: time_min=480, mem_mb=20000, cpus=16
    wrapper:
        f"{wrappers_version}/bio/sambamba/sort"

rule sambamba_sort_pe:
    input:
        "hisat2/{trimmer}_pe/{sample}.{unit}.bam"
    output:
        "hisat2/{trimmer}_pe/{sample}.{unit}.sorted.bam"
    log:
        "logs/hisat2/{trimmer}_pe/sambamba-sort/{sample}.{unit}.log"
    params: ""
    threads: 16 
    wildcard_constraints:
        unit=r"rep\d+"
    resources: time_min=480, mem_mb=20000, cpus=16
    wrapper:
        f"{wrappers_version}/bio/sambamba/sort"

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
        unit=r"rep\d+"
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

##Salmon Section
rule salmon_index:
    input:
        sequences=config["ref"]["transcriptomefa"],
    output:
        multiext(
            "salmon/transcriptome_index/",
            "complete_ref_lens.bin",
            "ctable.bin",
            "ctg_offsets.bin",
            "duplicate_clusters.tsv",
            "info.json",
            "mphf.bin",
            "pos.bin",
            "pre_indexing.log",
            "rank.bin",
            "refAccumLengths.bin",
            "ref_indexing.log",
            "reflengths.bin",
            "refseq.bin",
            "seq.bin",
            "versionInfo.json",
        ),
        idxdir=directory("salmon/transcriptome_index"),
    log:
        "logs/salmon/transcriptome_index.log",
    threads: 16
    resources: time_min=320, mem_mb=40000, cpus=16
    params:
        # optional parameters
        extra="",
    wrapper:
        f"{wrappers_version}/bio/salmon/index"


rule salmon_quant_reads_pe:
    input:
        r1="trimmed/{trimmer}_pe/{sample}.{unit}.1.fastq.gz",
        r2="trimmed/{trimmer}_pe/{sample}.{unit}.2.fastq.gz",
        index="salmon/transcriptome_index",
    output:
        quant="salmon/{trimmer}_pe/{sample}.{unit}/quant.sf",
        lib="salmon/{trimmer}_pe/{sample}.{unit}/lib_format_counts.json",
        dirout=directory("salmon/{trimmer}_pe/{sample}.{unit}/"),
    log:
        "logs/salmon/{trimmer}_pe/{sample}.{unit}.log",
    params:
        # optional parameters
        libtype="A",
        extra="--gcBias",
    threads: 16
    resources: time_min=320, mem_mb=40000, cpus=16
    wrapper:
        f"{wrappers_version}/bio/salmon/quant"

rule salmon_quant_reads_se:
    input:
        r="trimmed/{trimmer}_se/{sample}.{unit}.1.fastq.gz",
        index="salmon/transcriptome_index/",
    output:
        quant="salmon/{trimmer}_se/{sample}.{unit}/quant.sf",
        lib="salmon/{trimmer}_se/{sample}.{unit}/lib_format_counts.json",
        dirout=directory("salmon/{trimmer}_se/{sample}.{unit}/"),
    log:
        "logs/salmon/{trimmer}_se/{sample}.{unit}.log",
    params:
        # optional parameters
        libtype="A",
        extra="--gcBias",
    threads: 16
    resources: time_min=320, mem_mb=40000, cpus=16
    wrapper:
        f"{wrappers_version}/bio/salmon/quant"

#rule symlink_salmon_quants:
#    input:
#        quant="salmon/{trimmer}_{pese}/{sample}.{unit}/quant.sf",
#        lib="salmon/{trimmer}_{pese}/{sample}.{unit}/lib_format_counts.json",
#    output:
#        quant="salmon/{trimmer}_{pese}/{sample}.{unit}/{sample}.{unit}.quant.sf",
#        lib="salmon/{trimmer}_{pese}/{sample}.{unit}/{sample}.{unit}.lib_format_counts.json",
#    params:
#        #quant="{sample}.{unit}.quant.sf",
#        #lib="{sample}.{unit}.lib_format_counts.json",
#        quant="../../../salmon/{trimmer}_{pese}/{sample}.{unit}/quant.sf",
#        lib="../../../salmon/{trimmer}_{pese}/{sample}.{unit}/lib_format_counts.json",
#    threads: 1
#    resources: time_min=320, mem_mb=2000, cpus=1
#    shell:
#        "ln -s {params.quant} {output.quant} && ln -s {params.lib} {output.lib}"


rule salmon_merge_quants:
    input:
        quant=expand("salmon/{trimmer}_{pese}/{unit.sample}.{unit.unit}/", trimmer=trimmers, pese=pese, unit=units.itertuples())
    output:
        "salmon/{trimmer}_{pese}/merged_quant.tsv",
    threads: 1
    resources: time_min=320, mem_mb=20000, cpus=1
    conda: "../envs/salmon.yaml"
    log:
        "logs/salmon/{trimmer}_{pese}/merge_quant.log",
    shell:
        "salmon quantmerge --quants {input.quant} -o {output} 2> {log}"

