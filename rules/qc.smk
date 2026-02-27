def get_fastq(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

def get_fastq1(wildcards):
    fq1 = units.loc[(wildcards.sample, wildcards.unit), ["fq1"]].dropna()
    return fq1

def get_fastq2(wildcards):
    fq2 = units.loc[(wildcards.sample, wildcards.unit), ["fq2"]].dropna()
    return fq2

def get_picard_input(wildcards):
    if wildcards.aligner == "star":
        return f"star/{wildcards.trimmer}_{wildcards.pese}/{wildcards.sample}.{wildcards.unit}/{wildcards.sample}.{wildcards.unit}_Aligned.sortedByCoord.out.bam"
    elif wildcards.aligner == "hisat2":
        return f"hisat2/{wildcards.trimmer}_{wildcards.pese}/{wildcards.sample}.{wildcards.unit}.sorted.bam"

rule picard_mark_duplicates:
    input:
        bams = get_picard_input
    output:
        bam = "{aligner}/{trimmer}_{pese}/{sample}.{unit}.marked.bam",
        metrics = "{aligner}/{trimmer}_{pese}/{sample}.{unit}.marked.metrics.txt"
    log:
        "logs/picard/mark_duplicates/{sample}.{unit}.{aligner}.{trimmer}_{pese}.log"
    params:
        extra = "--REMOVE_DUPLICATES false"
    resources: 
        time_min=320, 
        mem_mb=40000, 
        cpus=1
    wrapper:
        f"{wrappers_version}/bio/picard/markduplicates"


#rule gtf_to_refflat:
#    input:
#        gtf=config["ref"]["annotation"] 
#    output:
#        # Uses f-string to dynamically insert the index name
#        refflat=f"refs/{config['ref']['index']}.refFlat"
#    log:
#        # It's also good practice to make the log file dynamic
#        f"logs/refs/gtfToGenePred_{config['ref']['index']}.log"
#    params:
#        convert_out="PicardCollectRnaSeqMetrics",
#        extra="-genePredExt -geneNameAsName2"
#    wrapper:
#        f"{wrappers_version}/bio/ucsc/gtfToGenePred"
rule gtf_to_refflat:
    input:
        gtf=config["ref"]["annotation"]
    output:
        refflat=f"refs/{config['ref']['index']}.refFlat"
    log:
        f"logs/refs/gtfToGenePred_{config['ref']['index']}.log"
    conda:
        "../envs/ucsc.yaml" 
    shell:
        """
        # 1. Create a standard 10-column genePred safely
        gtfToGenePred -ignoreGroupsWithoutExons {input.gtf} {output.refflat}.tmp 2> {log}
        
        # 2. Duplicate the first column to create the 11-column refFlat Picard demands
        awk 'BEGIN {{OFS="\\t"}} {{print $1, $0}}' {output.refflat}.tmp > {output.refflat}
        
        # 3. Clean up
        rm {output.refflat}.tmp
        """

rule picard_create_dict:
    input:
        config["ref"]["genomefa"]
    output:
        f"refs/{config['ref']['index']}.dict"
    log:
        f"logs/refs/create_dict_{config['ref']['index']}.log"
    params:
        extra=""
    resources:
        mem_mb=4000
    wrapper:
        f"{wrappers_version}/bio/picard/createsequencedictionary"

rule find_rrna_barrnap:
    input:
        fasta=config["ref"]["genomefa"]
    output:
        gff=f"refs/{config['ref']['index']}.barrnap.gff",
        fasta=f"refs/{config['ref']['index']}.barrnap_hits.fasta"
    params:
        kingdom="euk",
        extra=""
    resources: 
        time_min=320, 
        mem_mb=40000, 
        cpus=8
    log:
        f"logs/refs/barrnap_{config['ref']['index']}.log"
    wrapper:
        f"{wrappers_version}/bio/barrnap"

rule barrnap_gff_to_bed:
    input:
        f"refs/{config['ref']['index']}.barrnap.gff"
    output:
        f"refs/{config['ref']['index']}.barrnap.bed"
    shell:
        "grep -v '^#' {input} | awk -v OFS='\\t' '{{print $1, $4-1, $5, $9, \"0\", $7}}' > {output}"

rule picard_bed_to_interval_list:
    input:
        bed=f"refs/{config['ref']['index']}.barrnap.bed",
        dict=f"refs/{config['ref']['index']}.dict"
    output:
        f"refs/{config['ref']['index']}.rRNA.interval_list"
    log:
        f"logs/refs/bedtointervallist_{config['ref']['index']}.log"
    params:
        extra="--SORT true"
    resources: 
        time_min=320, 
        mem_mb=4000, 
        cpus=1
    wrapper:
        f"{wrappers_version}/bio/picard/bedtointervallist"

rule filter_discordant_for_picard:
    input:
        bam="{aligner}/{trimmer}_{pese}/{sample}.{unit}.marked.bam"
    output:
        bam=temp("qc/picard/{aligner}/{trimmer}_{pese}/{sample}.{unit}.concordant.bam")
    log:
        "logs/picard/filter_discordant/{sample}.{unit}.{aligner}.{trimmer}_{pese}.log"
    resources: 
        time_min=320, 
        mem_mb=4000, 
        cpus=4
    conda:
        "../envs/hisat2.yaml" 
    threads: 4
    shell:
        "samtools view -h -@ {threads} {input.bam} | awk '$1 ~ /^@/ || $7 == \"=\" || $7 == \"*\"' | samtools view -b -@ {threads} - > {output.bam} 2> {log}"
         
rule picard_rnaseq_metrics:
    input:
        bam="qc/picard/{aligner}/{trimmer}_{pese}/{sample}.{unit}.concordant.bam",
        refflat=f"refs/{config['ref']['index']}.refFlat",
        ribo=f"refs/{config['ref']['index']}.rRNA.interval_list"
    output:
        "qc/picard/{aligner}/{trimmer}_{pese}/{sample}.{unit}.rnaseq_metrics.txt"
    log:
        "logs/picard/rnaseq_metrics/{sample}.{unit}.{aligner}.{trimmer}_{pese}.log"
    params:
        strand="NONE", 
        extra=lambda wildcards, input: f"--VALIDATION_STRINGENCY SILENT --RIBOSOMAL_INTERVALS {input.ribo}"
    resources:
        time_min=120,
        mem_mb=8000, 
        cpus=1
    wrapper: 
        f"{wrappers_version}/bio/picard/collectrnaseqmetrics"

rule rseqc_stat_star:
    input:
        "star/{trimmer}/{sample}.{unit}/{sample}.{unit}_Aligned.sortedByCoord.out.bam",
        "star/{trimmer}/{sample}.{unit}/{sample}.{unit}_Aligned.sortedByCoord.out.bam.bai",
    output:
        "qc/star/{trimmer}/rseqc/{sample}.{unit}.stats.txt",
    log:
        "logs/star/{trimmer}/rseqc/rseqc_stat/{sample}.{unit}.log",
    conda:
        "../envs/rseqc.yaml"
    resources: time_min=320, mem_mb=20000, cpus=1
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"

rule rseqc_stat_hisat2:
    input:
        "hisat2/{trimmer}/{sample}.{unit}.sorted.bam",
        "hisat2/{trimmer}/{sample}.{unit}.sorted.bam.bai",
    output:
        "qc/hisat2/{trimmer}/rseqc/{sample}.{unit}.stats.txt",
    log:
        "logs/hisat2/{trimmer}/rseqc/rseqc_stat/{sample}.{unit}.log",
    conda:
        "../envs/rseqc.yaml"
    resources: time_min=320, mem_mb=20000, cpus=1
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"




rule rseqc_readdup_star:
    input:
        "star/{trimmer}/{sample}.{unit}/{sample}.{unit}_Aligned.sortedByCoord.out.bam",
        "star/{trimmer}/{sample}.{unit}/{sample}.{unit}_Aligned.sortedByCoord.out.bam.bai",
    output:
        "qc/star/{trimmer}/rseqc/{sample}.{unit}.readdup.DupRate_plot.pdf",
    log:
        "logs/star/{trimmer}/rseqc/rseqc_readdup/{sample}.{unit}.log",
    params:
        prefix=lambda w, output: strip_suffix(output[0], ".DupRate_plot.pdf"),
    resources: time_min=320, mem_mb=40000, cpus=2
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_duplication.py -i {input} -o {params.prefix} > {log} 2>&1"

rule rseqc_readdup_hisat2:
    input:
        "hisat2/{trimmer}/{sample}.{unit}.sorted.bam",
        "hisat2/{trimmer}/{sample}.{unit}.sorted.bam.bai",
    output:
        "qc/hisat2/{trimmer}/rseqc/{sample}.{unit}.readdup.DupRate_plot.pdf",
    log:
        "logs/hisat2/{trimmer}/rseqc/rseqc_readdup/{sample}.{unit}.log",
    params:
        prefix=lambda w, output: strip_suffix(output[0], ".DupRate_plot.pdf"),
    resources: time_min=320, mem_mb=40000, cpus=2
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_duplication.py -i {input} -o {params.prefix} > {log} 2>&1"


rule rseqc_readgc_star:
    input:
        "star/{trimmer}/{sample}.{unit}/{sample}.{unit}_Aligned.sortedByCoord.out.bam",
        "star/{trimmer}/{sample}.{unit}/{sample}.{unit}_Aligned.sortedByCoord.out.bam.bai",
    output:
        "qc/star/{trimmer}/rseqc/{sample}.{unit}.readgc.GC_plot.pdf",
    log:
        "logs/star/{trimmer}/rseqc/rseqc_readgc/{sample}.{unit}.log",
    resources: time_min=320, mem_mb=20000, cpus=1
    params:
        prefix=lambda w, output: strip_suffix(output[0], ".GC_plot.pdf"),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_GC.py -i {input} -o {params.prefix} > {log} 2>&1"

rule rseqc_readgc_hisat2:
    input:
        "hisat2/{trimmer}/{sample}.{unit}.sorted.bam",
        "hisat2/{trimmer}/{sample}.{unit}.sorted.bam.bai",
    output:
        "qc/hisat2/{trimmer}/rseqc/{sample}.{unit}.readgc.GC_plot.pdf",
    log:
        "logs/hisat2/{trimmer}/rseqc/rseqc_readgc/{sample}.{unit}.log",
    resources: time_min=320, mem_mb=20000, cpus=1
    params:
        prefix=lambda w, output: strip_suffix(output[0], ".GC_plot.pdf"),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_GC.py -i {input} -o {params.prefix} > {log} 2>&1"

rule fastqc_pretrim_r1:
    input:
       get_fastq1
    output:
        html="qc/fastqc_pretrim/{sample}.{unit}_r1.html",
        zip="qc/fastqc_pretrim/{sample}.{unit}_r1_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc_pretrim/{sample}.{unit}_r1.log"
    resources: time_min=320, mem_mb=10000, cpus=1
    threads: 1
    wrapper:
        f"{wrappers_version}/bio/fastqc"

rule fastqc_pretrim_r2:
    input:
       get_fastq2
    output:
        html="qc/fastqc_pretrim/{sample}.{unit}_r2.html",
        zip="qc/fastqc_pretrim/{sample}.{unit}_r2_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc_pretrim/{sample}.{unit}_r2.log"
    resources: time_min=320, mem_mb=10000, cpus=1
    threads: 1
    wrapper:
        f"{wrappers_version}/bio/fastqc"

rule fastqc_posttrim_r1:
    input:
         "trimmed/{trimmer}_{pese}/{sample}.{unit}.1.fastq.gz"
    output:
        html="qc/fastqc_posttrim/{trimmer}_{pese}/{sample}.{unit}_r1.html",
        zip="qc/fastqc_posttrim/{trimmer}_{pese}/{sample}.{unit}_r1_fastqc.zip" 
    params: ""
    log:
        "logs/fastqc_posttrim/{trimmer}_{pese}/{sample}.{unit}_r1.log"
    resources: time_min=320, mem_mb=10000, cpus=1
    threads: 1
    wrapper:
        f"{wrappers_version}/bio/fastqc"

rule fastqc_posttrim_r2:
    input:
        "trimmed/{trimmer}_{pese}/{sample}.{unit}.2.fastq.gz"
    output:
        html="qc/fastqc_posttrim/{trimmer}_{pese}/{sample}.{unit}_r2.html",
        zip="qc/fastqc_posttrim/{trimmer}_{pese}/{sample}.{unit}_r2_fastqc.zip"
    params: ""
    log:
        "logs/fastqc_posttrim/{trimmer}_{pese}/{sample}.{unit}_r2.log"
    resources: time_min=320, mem_mb=10000, cpus=1
    threads: 1
    wrapper:
        f"{wrappers_version}/bio/fastqc"

#rule fastqc_posttrim_r2:
#    input:
#        "trimmed/{trimmer}/{sample}.{unit}.2.fastq.gz"
#    output:
#        html="qc/fastqc_posttrim/{trimmer}/{sample}.{unit}_r2.html",
#        zip="qc/fastqc_posttrim/{trimmer}/{sample}.{unit}_r2_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
#    params: ""
#    log:
#        "logs/fastqc_posttrim/{trimmer}/{sample}.{unit}_r2.log"
#    resources: time_min=320, mem_mb=20000, cpus=1
#    threads: 1
#    wrapper:
#        #"v0.75.0/bio/fastqc"
#        f"{wrappers_version}/bio/fastqc"

rule multiqc_pre_pe:
    input:
        expand("qc/fastqc_pretrim/{unit.sample}.{unit.unit}_r1_fastqc.zip", unit=units.itertuples(), trimmer=trimmers),
        expand("qc/fastqc_pretrim/{unit.sample}.{unit.unit}_r2_fastqc.zip", unit=units.itertuples(), trimmer=trimmers)
    output:
        "qc/multiqc_report_pretrim_pe.html"
    log:
        "logs/multiqc_pretrim_pe.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        f"{wrappers_version}/bio/multiqc"

rule multiqc_pre_se:
    input:
        expand("qc/fastqc_pretrim/{unit.sample}.{unit.unit}_r1_fastqc.zip", unit=units.itertuples(), trimmer=trimmers),
    output:
        "qc/multiqc_report_pretrim_se.html"
    log:
        "logs/multiqc_pretrim_se.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        f"{wrappers_version}/bio/multiqc"

rule multiqc_post_trimmomatic:
    input:
        expand("logs/trimmomatic/{unit.sample}.{unit.unit}.log", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/trimmomatic/{unit.sample}.{unit.unit}_r1_fastqc.zip", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/trimmomatic/{unit.sample}.{unit.unit}_r2_fastqc.zip", unit=units.itertuples())
    output:
        "qc/multiqc_report_posttrim_trimmomatic.html"
    log:
        "logs/multiqc_posttrim_trimmomatic.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        f"{wrappers_version}/bio/multiqc"

rule multiqc_post_fastp:
    input:
        expand("report/fastp/{unit.sample}.{unit.unit}.fastp.json", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/fastp/{unit.sample}.{unit.unit}_r1_fastqc.zip", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/fastp/{unit.sample}.{unit.unit}_r2_fastqc.zip", unit=units.itertuples())
    output:
        "qc/multiqc_report_posttrim_fastp.html"
    log:
        "logs/multiqc_posttrim_fastp.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        f"{wrappers_version}/bio/multiqc"

rule multiqc_post_trimgalore:
    input:
        expand("logs/trimgalore/{unit.sample}.{unit.unit}.log", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/trimgalore/{unit.sample}.{unit.unit}_r1_fastqc.zip", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/trimgalore/{unit.sample}.{unit.unit}_r2_fastqc.zip", unit=units.itertuples())
    output:
        "qc/multiqc_report_posttrim_trimgalore.html"
    log:
        "logs/multiqc_posttrim_trimgalore.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        f"{wrappers_version}/bio/multiqc"


rule multiqc_star_fastp_pe:
    input:
        expand("star/fastp_pe/{unit.sample}.{unit.unit}/{unit.sample}.{unit.unit}_Aligned.sortedByCoord.out.bam", unit=units.itertuples()),
        expand("star/fastp_pe/{unit.sample}.{unit.unit}/{unit.sample}.{unit.unit}_Aligned.sortedByCoord.out.bam.flagstat", unit=units.itertuples()),
        expand("qc/picard/star/fastp_pe/{unit.sample}.{unit.unit}.rnaseq_metrics.txt", unit=units.itertuples()),
        "results/star/all.star.fastp_pe.featureCounts.summary",
        expand("report/fastp_pe/{unit.sample}.{unit.unit}.fastp.json", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/fastp_pe/{unit.sample}.{unit.unit}_r1_fastqc.zip", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/fastp_pe/{unit.sample}.{unit.unit}_r2_fastqc.zip", unit=units.itertuples()),
        expand("star/fastp_pe/{unit.sample}.{unit.unit}.marked.metrics.txt", unit=units.itertuples())
    output:
        "qc/multiqc_report_star_fastp_pe.html"
    log:
        "logs/multiqc_star_fastp_pe.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        f"{wrappers_version}/bio/multiqc"

rule multiqc_star_fastp_se:
    input:
        expand("star/fastp_se/{unit.sample}.{unit.unit}/{unit.sample}.{unit.unit}_Aligned.sortedByCoord.out.bam", unit=units.itertuples()),
        expand("star/fastp_se/{unit.sample}.{unit.unit}/{unit.sample}.{unit.unit}_Aligned.sortedByCoord.out.bam.flagstat", unit=units.itertuples()),
        expand("qc/picard/star/fastp_se/{unit.sample}.{unit.unit}.rnaseq_metrics.txt", unit=units.itertuples()),
        "results/star/all.star.fastp_se.featureCounts.summary",
        expand("report/fastp_se/{unit.sample}.{unit.unit}.fastp.json", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/fastp_se/{unit.sample}.{unit.unit}_r1_fastqc.zip", unit=units.itertuples()),
        expand("star/fastp_se/{unit.sample}.{unit.unit}.marked.metrics.txt", unit=units.itertuples())
    output:
        "qc/multiqc_report_star_fastp_se.html"
    log:
        "logs/multiqc_star_fastp_se.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        f"{wrappers_version}/bio/multiqc"

rule multiqc_star_trimgalore_pe:
    input:
        expand("star/trimgalore_pe/{unit.sample}.{unit.unit}/{unit.sample}.{unit.unit}_Aligned.sortedByCoord.out.bam", unit=units.itertuples()),
        expand("star/trimgalore_pe/{unit.sample}.{unit.unit}/{unit.sample}.{unit.unit}_Aligned.sortedByCoord.out.bam.flagstat", unit=units.itertuples()),
        "results/star/all.star.trimgalore_pe.featureCounts.summary",
        expand("trimmed/trimgalore_pe/{unit.sample}.{unit.unit}.1_trimming_report.txt", unit=units.itertuples()),
        expand("trimmed/trimgalore_pe/{unit.sample}.{unit.unit}.2_trimming_report.txt", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/trimgalore_pe/{unit.sample}.{unit.unit}_r1_fastqc.zip", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/trimgalore_pe/{unit.sample}.{unit.unit}_r2_fastqc.zip", unit=units.itertuples()),
        expand("star/trimgalore_pe/{unit.sample}.{unit.unit}.marked.metrics.txt", unit=units.itertuples()),
    output:
        "qc/multiqc_report_star_trimgalore_pe.html"
    log:
        "logs/multiqc_star_trimgalore_pe.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        f"{wrappers_version}/bio/multiqc"

rule multiqc_star_trimgalore_se:
    input:
        expand("star/trimgalore_se/{unit.sample}.{unit.unit}/{unit.sample}.{unit.unit}_Aligned.sortedByCoord.out.bam", unit=units.itertuples()),
        expand("star/trimgalore_se/{unit.sample}.{unit.unit}/{unit.sample}.{unit.unit}_Aligned.sortedByCoord.out.bam.flagstat", unit=units.itertuples()),
        "results/star/all.star.trimgalore_se.featureCounts.summary",
        expand("trimmed/trimgalore_se/{unit.sample}.{unit.unit}.1_trimming_report.txt", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/trimgalore_se/{unit.sample}.{unit.unit}_r1_fastqc.zip", unit=units.itertuples()),
        expand("star/trimgalore_se/{unit.sample}.{unit.unit}.marked.metrics.txt", unit=units.itertuples()),
    output:
        "qc/multiqc_report_star_trimgalore_se.html"
    log:
        "logs/multiqc_star_trimgalore_se.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        f"{wrappers_version}/bio/multiqc"


rule multiqc_hisat2_fastp_se:
    input:
        expand("hisat2/fastp_se/{unit.sample}.{unit.unit}.sorted.bam", unit=units.itertuples()),
        expand("hisat2/fastp_se/{unit.sample}.{unit.unit}.sorted.bam.flagstat", unit=units.itertuples()),
        expand("qc/picard/hisat2/fastp_se/{unit.sample}.{unit.unit}.rnaseq_metrics.txt", unit=units.itertuples()),
        "results/hisat2/all.hisat2.fastp_se.featureCounts.summary",
        #expand("qc/hisat2/fastp/rseqc/{unit.sample}.{unit.unit}.stats.txt", unit=units.itertuples()),
        #expand("qc/hisat2/fastp/rseqc/{unit.sample}.{unit.unit}.readdup.DupRate_plot.pdf", unit=units.itertuples()),
        #expand("qc/hisat2/fastp/rseqc/{unit.sample}.{unit.unit}.readgc.GC_plot.pdf", unit=units.itertuples()),
        expand("report/fastp_se/{unit.sample}.{unit.unit}.fastp.json", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/fastp_se/{unit.sample}.{unit.unit}_r1_fastqc.zip", unit=units.itertuples()),
        expand("hisat2/fastp_se/{unit.sample}.{unit.unit}.marked.metrics.txt", unit=units.itertuples()),
    output:
        "qc/multiqc_report_hisat2_fastp_se.html"
    log:
        "logs/multiqc_hisat2_fastp_se.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        f"{wrappers_version}/bio/multiqc"

rule multiqc_hisat2_fastp_pe:
    input:
        expand("hisat2/fastp_pe/{unit.sample}.{unit.unit}.sorted.bam", unit=units.itertuples()),
        expand("hisat2/fastp_pe/{unit.sample}.{unit.unit}.sorted.bam.flagstat", unit=units.itertuples()),
        expand("qc/picard/hisat2/fastp_pe/{unit.sample}.{unit.unit}.rnaseq_metrics.txt", unit=units.itertuples()),
        "results/hisat2/all.hisat2.fastp_pe.featureCounts.summary",
        expand("report/fastp_pe/{unit.sample}.{unit.unit}.fastp.json", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/fastp_pe/{unit.sample}.{unit.unit}_r1_fastqc.zip", unit=units.itertuples()),
        expand("hisat2/fastp_pe/{unit.sample}.{unit.unit}.marked.metrics.txt", unit=units.itertuples()),
    output:
        "qc/multiqc_report_hisat2_fastp_pe.html"
    log:
        "logs/multiqc_hisat2_fastp_pe.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        f"{wrappers_version}/bio/multiqc"


rule multiqc_hisat2_fastp_nofct:
    input:
        expand("hisat2/fastp_{pese}/{unit.sample}.{unit.unit}.sorted.bam", unit=units.itertuples(), pese=pese),
        expand("report/fastp_{pese}/{unit.sample}.{unit.unit}.fastp.json", unit=units.itertuples(), pese=pese),
        expand("hisat2/fastp_{pese}/{unit.sample}.{unit.unit}.sorted.bam.flagstat", unit=units.itertuples(), pese=pese),
        expand("hisat2/fastp_{pese}/{unit.sample}.{unit.unit}.marked.metrics.txt", unit=units.itertuples(), pese=pese),
        #expand("qc/fastqc_posttrim/fastp/{unit.sample}.{unit.unit}_r1_fastqc.zip", unit=units.itertuples()),
        #expand("qc/fastqc_posttrim/fastp/{unit.sample}.{unit.unit}_r2_fastqc.zip", unit=units.itertuples())
    output:
        "qc/multiqc_report_hisat2_fastp_nofct.html"
    log:
        "logs/multiqc_hisat2_fastp_nofct.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        f"{wrappers_version}/bio/multiqc"

rule multiqc_hisat2_trimgalore_pe:
    input:
        expand("hisat2/trimgalore_pe/{unit.sample}.{unit.unit}.sorted.bam", unit=units.itertuples()),
        expand("hisat2/trimgalore_pe/{unit.sample}.{unit.unit}.sorted.bam.flagstat", unit=units.itertuples()),
        "results/hisat2/all.hisat2.trimgalore_pe.featureCounts.summary",
        expand("trimmed/trimgalore_pe/{unit.sample}.{unit.unit}.1_trimming_report.txt", unit=units.itertuples()),
        expand("trimmed/trimgalore_pe/{unit.sample}.{unit.unit}.2_trimming_report.txt", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/trimgalore_pe/{unit.sample}.{unit.unit}_r1_fastqc.zip", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/trimgalore_pe/{unit.sample}.{unit.unit}_r2_fastqc.zip", unit=units.itertuples()),
        expand("hisat2/trimgalore_pe/{unit.sample}.{unit.unit}.marked.metrics.txt", unit=units.itertuples()),
    output:
        "qc/multiqc_report_hisat2_trimgalore_pe.html"
    log:
        "logs/multiqc_hisat2_trimgalore_pe.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        f"{wrappers_version}/bio/multiqc"

rule multiqc_hisat2_trimgalore_se:
    input:
        expand("hisat2/trimgalore_se/{unit.sample}.{unit.unit}.sorted.bam", unit=units.itertuples()),
        expand("hisat2/trimgalore_se/{unit.sample}.{unit.unit}.sorted.bam.flagstat", unit=units.itertuples()),
        "results/hisat2/all.hisat2.trimgalore_se.featureCounts.summary",
        expand("trimmed/trimgalore_se/{unit.sample}.{unit.unit}.1_trimming_report.txt", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/trimgalore_se/{unit.sample}.{unit.unit}_r1_fastqc.zip", unit=units.itertuples()),
        expand("hisat2/trimgalore_se/{unit.sample}.{unit.unit}.marked.metrics.txt", unit=units.itertuples()),
    output:
        "qc/multiqc_report_hisat2_trimgalore_se.html"
    log:
        "logs/multiqc_hisat2_trimgalore_se.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        f"{wrappers_version}/bio/multiqc"

rule multiqc_salmon_fastp_pe:
    input:
        expand("salmon/fastp_pe/{unit.sample}.{unit.unit}/quant.sf", unit=units.itertuples()),
        expand("report/fastp_pe/{unit.sample}.{unit.unit}.fastp.json", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/fastp_pe/{unit.sample}.{unit.unit}_r1_fastqc.zip", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/fastp_pe/{unit.sample}.{unit.unit}_r2_fastqc.zip", unit=units.itertuples())
    output:
        "qc/multiqc_report_salmon_fastp_pe.html"
    log:
        "logs/multiqc_salmon_fastp_pe.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        f"{wrappers_version}/bio/multiqc"

rule multiqc_salmon_trimgalore_pe:
    input:
        expand("salmon/trimgalore_pe/{unit.sample}.{unit.unit}/quant.sf", unit=units.itertuples()),
        expand("trimmed/trimgalore_pe/{unit.sample}.{unit.unit}.1_trimming_report.txt", unit=units.itertuples()),
        expand("trimmed/trimgalore_pe/{unit.sample}.{unit.unit}.2_trimming_report.txt", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/trimgalore_pe/{unit.sample}.{unit.unit}_r1_fastqc.zip", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/trimgalore_pe/{unit.sample}.{unit.unit}_r2_fastqc.zip", unit=units.itertuples())
    output:
        "qc/multiqc_report_salmon_trimgalore_pe.html"
    log:
        "logs/multiqc_salmon_trimgalore_pe.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        f"{wrappers_version}/bio/multiqc"

rule multiqc_salmon_fastp_se:
    input:
        expand("salmon/fastp_se/{unit.sample}.{unit.unit}/quant.sf", unit=units.itertuples()),
        expand("report/fastp_se/{unit.sample}.{unit.unit}.fastp.json", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/fastp_se/{unit.sample}.{unit.unit}_r1_fastqc.zip", unit=units.itertuples()),
    output:
        "qc/multiqc_report_salmon_fastp_se.html"
    log:
        "logs/multiqc_salmon_fastp_se.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        f"{wrappers_version}/bio/multiqc"

rule multiqc_salmon_trimgalore_se:
    input:
        expand("salmon/trimgalore_se/{unit.sample}.{unit.unit}/quant.sf", unit=units.itertuples()),
        expand("trimmed/trimgalore_se/{unit.sample}.{unit.unit}.1_trimming_report.txt", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/trimgalore_se/{unit.sample}.{unit.unit}_r1_fastqc.zip", unit=units.itertuples()),
    output:
        "qc/multiqc_report_salmon_trimgalore_se.html"
    log:
        "logs/multiqc_salmon_trimgalore_se.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        f"{wrappers_version}/bio/multiqc"


