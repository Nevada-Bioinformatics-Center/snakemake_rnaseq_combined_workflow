def get_fastq(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

def get_fastq1(wildcards):
    fq1 = units.loc[(wildcards.sample, wildcards.unit), ["fq1"]].dropna()
    return fq1

def get_fastq2(wildcards):
    fq2 = units.loc[(wildcards.sample, wildcards.unit), ["fq2"]].dropna()
    return fq2

## RSEQC


#rule rseqc_gtf2bed:
#    input:
#        config["ref"]["annotation"],
#    output:
#        bed="qc/rseqc/annotation.bed",
#        db=temp("qc/rseqc/annotation.db"),
#    log:
#        "logs/rseqc_gtf2bed.log",
#    resources: time_min=320, mem_mb=8000, cpus=1
#    conda:
#        "../envs/gffutils.yaml"
#    script:
#        "../scripts/gtf2bed.py"
#

#rule rseqc_junction_annotation:
#    input:
#        bam="star/{sample}-{unit}/{sample}.{unit}_Aligned.sortedByCoord.out.bam",
#        bai="star/{sample}-{unit}/{sample}.{unit}_Aligned.sortedByCoord.out.bam.bai",
#        bed="qc/rseqc/annotation.bed",
#    output:
#        "qc/rseqc/{sample}-{unit}.junctionanno.junction.bed",
#    priority: 1
#    log:
#        "logs/rseqc/rseqc_junction_annotation/{sample}-{unit}.log",
#    params:
#        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
#        prefix=lambda w, output: strip_suffix(output[0], ".junction.bed"),
#    resources: time_min=320, mem_mb=8000, cpus=1
#    conda:
#        "../envs/rseqc.yaml"
#    shell:
#        "junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
#        "> {log[0]} 2>&1"


#rule rseqc_junction_saturation:
#    input:
#        bam="star/{sample}-{unit}/{sample}.{unit}_Aligned.sortedByCoord.out.bam",
#        bai="star/{sample}-{unit}/{sample}.{unit}_Aligned.sortedByCoord.out.bam.bai",
#        bed="qc/rseqc/annotation.bed",
#    output:
#        "qc/rseqc/{sample}-{unit}.junctionsat.junctionSaturation_plot.pdf",
#    priority: 1
#    log:
#        "logs/rseqc/rseqc_junction_saturation/{sample}-{unit}.log",
#    params:
#        extra=r"-q 255",
#        prefix=lambda w, output: strip_suffix(
#            output[0], ".junctionSaturation_plot.pdf"
#        ),
#    resources: time_min=320, mem_mb=8000, cpus=1
#    conda:
#        "../envs/rseqc.yaml"
#    shell:
#        "junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
#        "> {log} 2>&1"
#


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


#rule rseqc_infer:
#    input:
#        bam="star/{sample}-{unit}/{sample}.{unit}_Aligned.sortedByCoord.out.bam",
#        bai="star/{sample}-{unit}/{sample}.{unit}_Aligned.sortedByCoord.out.bam.bai",
#        bed="qc/rseqc/annotation.bed",
#    output:
#        "qc/rseqc/{sample}-{unit}.infer_experiment.txt",
#    priority: 1
#    log:
#        "logs/rseqc/rseqc_infer/{sample}-{unit}.log",
#    resources: time_min=320, mem_mb=8000, cpus=1
#    conda:
#        "../envs/rseqc.yaml"
#    shell:
#        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"
#
#
#rule rseqc_innerdis:
#    input:
#        bam="star/{sample}-{unit}/{sample}.{unit}_Aligned.sortedByCoord.out.bam",
#        bai="star/{sample}-{unit}/{sample}.{unit}_Aligned.sortedByCoord.out.bam.bai",
#        bed="qc/rseqc/annotation.bed",
#    output:
#        "qc/rseqc/{sample}-{unit}.inner_distance_freq.inner_distance.txt",
#    priority: 1
#    log:
#        "logs/rseqc/rseqc_innerdis/{sample}-{unit}.log",
#    resources: time_min=320, mem_mb=8000, cpus=1
#    params:
#        prefix=lambda w, output: strip_suffix(output[0], ".inner_distance.txt"),
#    conda:
#        "../envs/rseqc.yaml"
#    shell:
#        "inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1"
#
#
#rule rseqc_readdis:
#    input:
#        bam="star/{sample}-{unit}/{sample}.{unit}_Aligned.sortedByCoord.out.bam",
#        bai="star/{sample}-{unit}/{sample}.{unit}_Aligned.sortedByCoord.out.bam.bai",
#        bed="qc/rseqc/annotation.bed",
#    output:
#        "qc/rseqc/{sample}-{unit}.readdistribution.txt",
#    priority: 1
#    log:
#        "logs/rseqc/rseqc_readdis/{sample}-{unit}.log",
#    resources: time_min=320, mem_mb=8000, cpus=1
#    conda:
#        "../envs/rseqc.yaml"
#    shell:
#        "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


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
        html="qc/fastqc_pretrim/{trimmer}/{sample}.{unit}_r1.html",
        zip="qc/fastqc_pretrim/{trimmer}/{sample}.{unit}_r1_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc_pretrim/{trimmer}/{sample}.{unit}_r1.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    threads: 1
    wrapper:
        "v0.75.0/bio/fastqc"

rule fastqc_pretrim_r2:
    input:
       get_fastq2
    output:
        html="qc/fastqc_pretrim/{trimmer}/{sample}.{unit}_r2.html",
        zip="qc/fastqc_pretrim/{trimmer}/{sample}.{unit}_r2_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc_pretrim/{trimmer}/{sample}.{unit}_r2.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    threads: 1
    wrapper:
        "v0.75.0/bio/fastqc"

rule fastqc_posttrim_r1:
    input:
        "trimmed/{trimmer}/{sample}.{unit}.1.fastq.gz"
    output:
        html="qc/fastqc_posttrim/{trimmer}/{sample}.{unit}_r1.html",
        zip="qc/fastqc_posttrim/{trimmer}/{sample}.{unit}_r1_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc_posttrim/{trimmer}/{sample}.{unit}_r1.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    threads: 1
    wrapper:
        "v0.75.0/bio/fastqc"

rule fastqc_posttrim_r2:
    input:
        "trimmed/{trimmer}/{sample}.{unit}.2.fastq.gz"
    output:
        html="qc/fastqc_posttrim/{trimmer}/{sample}.{unit}_r2.html",
        zip="qc/fastqc_posttrim/{trimmer}/{sample}.{unit}_r2_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc_posttrim/{trimmer}/{sample}.{unit}_r2.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    threads: 1
    wrapper:
        "v0.75.0/bio/fastqc"

rule multiqc_pre:
    input:
        expand("qc/fastqc_pretrim/{trimmer}/{unit.sample}.{unit.unit}_r1_fastqc.zip", unit=units.itertuples(), trimmer=trimmers),
        expand("qc/fastqc_pretrim/{trimmer}/{unit.sample}.{unit.unit}_r2_fastqc.zip", unit=units.itertuples(), trimmer=trimmers)
    output:
        "qc/multiqc_report_pretrim.html"
    log:
        "logs/multiqc_pre.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        #"0.84.0/bio/multiqc"
        f"{WRAPPER_PREFIX}/master/bio/multiqc"

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
        #"0.84.0/bio/multiqc"
        f"{WRAPPER_PREFIX}/master/bio/multiqc"

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
        #"0.84.0/bio/multiqc"
        f"{WRAPPER_PREFIX}/master/bio/multiqc"

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
        #"0.84.0/bio/multiqc"
        f"{WRAPPER_PREFIX}/master/bio/multiqc"

rule multiqc_star_trimmomatic:
    input:
        expand("star/trimmomatic/{unit.sample}.{unit.unit}/{unit.sample}.{unit.unit}_Aligned.sortedByCoord.out.bam", unit=units.itertuples()),
        #expand("results/featureCounts/{unit.sample}-{unit.unit}.featureCounts.summary", unit=units.itertuples()),
        "results/star/all.star.trimmomatic.featureCounts.summary",
        #expand("stats/{unit.sample}-{unit.unit}.isize.txt", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.junctionanno.junction.bed", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.junctionsat.junctionSaturation_plot.pdf", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.infer_experiment.txt", unit=units.itertuples()),
        expand("qc/star/trimmomatic/rseqc/{unit.sample}.{unit.unit}.stats.txt", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.inner_distance_freq.inner_distance.txt", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.readdistribution.txt", unit=units.itertuples()),
        expand("qc/star/trimmomatic/rseqc/{unit.sample}.{unit.unit}.readdup.DupRate_plot.pdf", unit=units.itertuples()),
        expand("qc/star/trimmomatic/rseqc/{unit.sample}.{unit.unit}.readgc.GC_plot.pdf", unit=units.itertuples()),
        #expand("logs/rseqc/rseqc_junction_annotation/{unit.sample}-{unit.unit}.log", unit=units.itertuples()),
        expand("logs/trimmomatic/{unit.sample}.{unit.unit}.log", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/trimmomatic/{unit.sample}.{unit.unit}_r1_fastqc.zip", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/trimmomatic/{unit.sample}.{unit.unit}_r2_fastqc.zip", unit=units.itertuples())
    output:
        "qc/multiqc_report_star_trimmomatic.html"
    log:
        "logs/multiqc_star_trimmomatic.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        #"0.84.0/bio/multiqc"
        f"{WRAPPER_PREFIX}/master/bio/multiqc"

rule multiqc_star_fastp:
    input:
        expand("star/fastp/{unit.sample}.{unit.unit}/{unit.sample}.{unit.unit}_Aligned.sortedByCoord.out.bam", unit=units.itertuples()),
        #expand("results/featureCounts/{unit.sample}-{unit.unit}.featureCounts.summary", unit=units.itertuples()),
        "results/star/all.star.fastp.featureCounts.summary",
        #expand("stats/{unit.sample}-{unit.unit}.isize.txt", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.junctionanno.junction.bed", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.junctionsat.junctionSaturation_plot.pdf", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.infer_experiment.txt", unit=units.itertuples()),
        expand("qc/star/fastp/rseqc/{unit.sample}.{unit.unit}.stats.txt", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.inner_distance_freq.inner_distance.txt", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.readdistribution.txt", unit=units.itertuples()),
        expand("qc/star/fastp/rseqc/{unit.sample}.{unit.unit}.readdup.DupRate_plot.pdf", unit=units.itertuples()),
        expand("qc/star/fastp/rseqc/{unit.sample}.{unit.unit}.readgc.GC_plot.pdf", unit=units.itertuples()),
        #expand("logs/rseqc/rseqc_junction_annotation/{unit.sample}-{unit.unit}.log", unit=units.itertuples()),
        expand("report/fastp/{unit.sample}.{unit.unit}.fastp.json", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/fastp/{unit.sample}.{unit.unit}_r1_fastqc.zip", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/fastp/{unit.sample}.{unit.unit}_r2_fastqc.zip", unit=units.itertuples())
    output:
        "qc/multiqc_report_star_fastp.html"
    log:
        "logs/multiqc_star_fastp.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        #"0.84.0/bio/multiqc"
        f"{WRAPPER_PREFIX}/master/bio/multiqc"

rule multiqc_star_trimgalore:
    input:
        expand("star/trimgalore/{unit.sample}.{unit.unit}/{unit.sample}.{unit.unit}_Aligned.sortedByCoord.out.bam", unit=units.itertuples()),
        #expand("results/featureCounts/{unit.sample}-{unit.unit}.featureCounts.summary", unit=units.itertuples()),
        "results/star/all.star.trimgalore.featureCounts.summary",
        #expand("stats/{unit.sample}-{unit.unit}.isize.txt", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.junctionanno.junction.bed", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.junctionsat.junctionSaturation_plot.pdf", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.infer_experiment.txt", unit=units.itertuples()),
        expand("qc/star/trimgalore/rseqc/{unit.sample}.{unit.unit}.stats.txt", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.inner_distance_freq.inner_distance.txt", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.readdistribution.txt", unit=units.itertuples()),
        expand("qc/star/trimgalore/rseqc/{unit.sample}.{unit.unit}.readdup.DupRate_plot.pdf", unit=units.itertuples()),
        expand("qc/star/trimgalore/rseqc/{unit.sample}.{unit.unit}.readgc.GC_plot.pdf", unit=units.itertuples()),
        #expand("logs/rseqc/rseqc_junction_annotation/{unit.sample}-{unit.unit}.log", unit=units.itertuples()),
        expand("logs/trimgalore/{unit.sample}.{unit.unit}.log", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/trimgalore/{unit.sample}.{unit.unit}_r1_fastqc.zip", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/trimgalore/{unit.sample}.{unit.unit}_r2_fastqc.zip", unit=units.itertuples())
    output:
        "qc/multiqc_report_star_trimgalore.html"
    log:
        "logs/multiqc_star_trimgalore.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        f"{WRAPPER_PREFIX}/master/bio/multiqc"
        #"0.84.0/bio/multiqc"

rule multiqc_hisat2_trimmomatic:
    input:
        expand("hisat2/trimmomatic/{unit.sample}.{unit.unit}.sorted.bam", unit=units.itertuples()),
        #expand("results/featureCounts/{unit.sample}-{unit.unit}.featureCounts.summary", unit=units.itertuples()),
        "results/hisat2/all.hisat2.trimmomatic.featureCounts.summary",
        #expand("stats/{unit.sample}-{unit.unit}.isize.txt", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.junctionanno.junction.bed", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.junctionsat.junctionSaturation_plot.pdf", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.infer_experiment.txt", unit=units.itertuples()),
        expand("qc/hisat2/trimmomatic/rseqc/{unit.sample}.{unit.unit}.stats.txt", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.inner_distance_freq.inner_distance.txt", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.readdistribution.txt", unit=units.itertuples()),
        expand("qc/hisat2/trimmomatic/rseqc/{unit.sample}.{unit.unit}.readdup.DupRate_plot.pdf", unit=units.itertuples()),
        expand("qc/hisat2/trimmomatic/rseqc/{unit.sample}.{unit.unit}.readgc.GC_plot.pdf", unit=units.itertuples()),
        #expand("logs/rseqc/rseqc_junction_annotation/{unit.sample}-{unit.unit}.log", unit=units.itertuples()),
        expand("logs/trimmomatic/{unit.sample}.{unit.unit}.log", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/trimmomatic/{unit.sample}.{unit.unit}_r1_fastqc.zip", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/trimmomatic/{unit.sample}.{unit.unit}_r2_fastqc.zip", unit=units.itertuples())
    output:
        "qc/multiqc_report_hisat2_trimmomatic.html"
    log:
        "logs/multiqc_hisat2_trimmomatic.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        #"0.84.0/bio/multiqc"
        f"{WRAPPER_PREFIX}/master/bio/multiqc"

rule multiqc_hisat2_fastp:
    input:
        expand("hisat2/fastp/{unit.sample}.{unit.unit}.sorted.bam", unit=units.itertuples()),
        #expand("results/featureCounts/{unit.sample}-{unit.unit}.featureCounts.summary", unit=units.itertuples()),
        "results/hisat2/all.hisat2.fastp.featureCounts.summary",
        #expand("stats/{unit.sample}-{unit.unit}.isize.txt", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.junctionanno.junction.bed", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.junctionsat.junctionSaturation_plot.pdf", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.infer_experiment.txt", unit=units.itertuples()),
        expand("qc/hisat2/fastp/rseqc/{unit.sample}.{unit.unit}.stats.txt", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.inner_distance_freq.inner_distance.txt", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.readdistribution.txt", unit=units.itertuples()),
        expand("qc/hisat2/fastp/rseqc/{unit.sample}.{unit.unit}.readdup.DupRate_plot.pdf", unit=units.itertuples()),
        expand("qc/hisat2/fastp/rseqc/{unit.sample}.{unit.unit}.readgc.GC_plot.pdf", unit=units.itertuples()),
        #expand("logs/rseqc/rseqc_junction_annotation/{unit.sample}-{unit.unit}.log", unit=units.itertuples()),
        expand("report/fastp/{unit.sample}.{unit.unit}.fastp.json", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/fastp/{unit.sample}.{unit.unit}_r1_fastqc.zip", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/fastp/{unit.sample}.{unit.unit}_r2_fastqc.zip", unit=units.itertuples())
    output:
        "qc/multiqc_report_hisat2_fastp.html"
    log:
        "logs/multiqc_hisat2_fastp.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        #"0.84.0/bio/multiqc"
        f"{WRAPPER_PREFIX}/master/bio/multiqc"

rule multiqc_hisat2_trimgalore:
    input:
        expand("hisat2/trimgalore/{unit.sample}.{unit.unit}.sorted.bam", unit=units.itertuples()),
        #expand("results/featureCounts/{unit.sample}-{unit.unit}.featureCounts.summary", unit=units.itertuples()),
        "results/hisat2/all.hisat2.trimgalore.featureCounts.summary",
        #expand("stats/{unit.sample}-{unit.unit}.isize.txt", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.junctionanno.junction.bed", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.junctionsat.junctionSaturation_plot.pdf", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.infer_experiment.txt", unit=units.itertuples()),
        expand("qc/hisat2/trimgalore/rseqc/{unit.sample}.{unit.unit}.stats.txt", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.inner_distance_freq.inner_distance.txt", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.readdistribution.txt", unit=units.itertuples()),
        expand("qc/hisat2/trimgalore/rseqc/{unit.sample}.{unit.unit}.readdup.DupRate_plot.pdf", unit=units.itertuples()),
        expand("qc/hisat2/trimgalore/rseqc/{unit.sample}.{unit.unit}.readgc.GC_plot.pdf", unit=units.itertuples()),
        #expand("logs/rseqc/rseqc_junction_annotation/{unit.sample}-{unit.unit}.log", unit=units.itertuples()),
        expand("logs/trimgalore/{unit.sample}.{unit.unit}.log", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/trimgalore/{unit.sample}.{unit.unit}_r1_fastqc.zip", unit=units.itertuples()),
        expand("qc/fastqc_posttrim/trimgalore/{unit.sample}.{unit.unit}_r2_fastqc.zip", unit=units.itertuples())
    output:
        "qc/multiqc_report_hisat2_trimgalore.html"
    log:
        "logs/multiqc_hisat2_trimgalore.log"
    resources: time_min=320, mem_mb=20000, cpus=1
    wrapper:
        #"0.84.0/bio/multiqc"
        f"{WRAPPER_PREFIX}/master/bio/multiqc"
