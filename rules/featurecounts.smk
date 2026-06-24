rule fix_gtf:
    input:
        gtf=config["ref"]["annotation"]
    output:
        gtf="results/ref/prepared_annotation.gtf"
    run:
        # Step 1: Scan the file to see if a fix is even necessary
        needs_fix = False
        with open(input.gtf, 'r') as f:
            for line in f:
                if 'gene_id ""' in line:
                    needs_fix = True
                    break
        
        # Step 2: Fix the file OR symlink it
        if needs_fix:
            print(f"Empty gene_ids found in {input.gtf}. Creating fixed version...")
            with open(input.gtf, 'r') as infile, open(output.gtf, 'w') as outfile:
                for line in infile:
                    if line.startswith("#"):
                        outfile.write(line)
                        continue
                    if 'gene_id ""' in line:
                        match = re.search(r'transcript_id "([^"]+)"', line)
                        if match:
                            line = line.replace('gene_id ""', f'gene_id "{match.group(1)}"')
                    outfile.write(line)
        else:
            print(f"GTF is clean. Symlinking {input.gtf}...")
            # Use absolute path to ensure the symlink doesn't break
            abs_input = os.path.abspath(input.gtf)
            os.symlink(abs_input, output.gtf)

# STAR featureCounts rule with multiext syntax
rule featurecounts_onefile_star:
    input:
        samples=lambda wc: expand(
            "star/{trimmer}_{pese}/{unit.sample}.{unit.unit}/{unit.sample}.{unit.unit}_Aligned.sortedByCoord.out.bam",
            unit=units.itertuples(), trimmer=wc.trimmer, pese=wc.pese
        ),
        annotation="results/ref/prepared_annotation.gtf",
       # annotation=config["ref"]["annotation"],
        fasta=config["ref"]["genomefa"]
    output:
        multiext(
            "results/star/all.star.{trimmer}_{pese}{fct_mode}",
            ".featureCounts",
            ".featureCounts.summary"
        )
    log:
        "logs/star/{trimmer}/featurecount_all_{pese}{fct_mode}.log"
    params:
        tmp_dir="",
        r_path="",
        extra=lambda wc: (
            ("-p -B -C " if wc.pese == "pe" else "") +
            ("-M --fraction " if wc.fct_mode == "_multifrac"
             else "-M "          if wc.fct_mode == "_multi"
             else "") +
            config["params"]["featurecounts"]
        ),
    threads: config["params"]["featurecountscpu"]
    resources:
        time_min=220,
        mem_mb=config["params"]["featurecountsram"],
        cpus=config["params"]["featurecountscpu"]
    wrapper:
        f"{wrappers_version}/bio/subread/featurecounts"

# HISAT2 featureCounts rule with multiext syntax
rule featurecounts_onefile_hisat2:
    input:
        samples=lambda wc: expand(
            "hisat2/{trimmer}_{pese}/{unit.sample}.{unit.unit}.sorted.bam",
            unit=units.itertuples(), trimmer=wc.trimmer, pese=wc.pese
        ),
        annotation="results/ref/prepared_annotation.gtf",
        #annotation=config["ref"]["annotation"],
        fasta=config["ref"]["genomefa"]
    output:
        multiext(
            "results/hisat2/all.hisat2.{trimmer}_{pese}{fct_mode}",
            ".featureCounts",
            ".featureCounts.summary"
        )
    log:
        "logs/hisat2/{trimmer}/featurecount_all_{pese}{fct_mode}.log"
    params:
        tmp_dir="",
        r_path="",
        extra=lambda wc: (
            ("-p -B -C " if wc.pese == "pe" else "") +
            ("-M --fraction " if wc.fct_mode == "_multifrac"
             else "-M "          if wc.fct_mode == "_multi"
             else "") +
            config["params"]["featurecounts"]
        ),
    threads: config["params"]["featurecountscpu"]
    resources:
        time_min=220,
        mem_mb=config["params"]["featurecountsram"],
        cpus=config["params"]["featurecountscpu"]
    wrapper:
        f"{wrappers_version}/bio/subread/featurecounts"


rule fix_featurecounts_general:
    input:
        "results/{aligner}/all.{aligner}.{trimmer}_{pese}.featureCounts",
    output:
        cwd+"results/{aligner}/all.{aligner}.{trimmer}_{pese}.fixcol2.featureCounts",
    log:
        "logs/star/fct_fix_{aligner}_{trimmer}_{pese}.log"
    threads: 1
    resources: time_min=220, mem_mb=8000, cpus=1
    shell: 
        "python3 scripts/fix_featurecounts_output.py -f {input} -c 2 > {output} 2> {log}"


rule fix_featurecounts_general_multi:
    input:
        "results/{aligner}/all.{aligner}.{trimmer}_{pese}_multi.featureCounts",
    output:
        cwd+"results/{aligner}/all.{aligner}.{trimmer}_{pese}_multi.fixcol2.featureCounts",
    log:
        "logs/star/fct_fix_{aligner}_{trimmer}_{pese}_multi.log"
    threads: 1
    resources: time_min=220, mem_mb=8000, cpus=1
    shell: 
        "python3 scripts/fix_featurecounts_output.py -f {input} -c 2 > {output} 2> {log}"

rule fix_featurecounts_general_multifrac:
    input:
        "results/{aligner}/all.{aligner}.{trimmer}_{pese}_multifrac.featureCounts",
    output:
        cwd+"results/{aligner}/all.{aligner}.{trimmer}_{pese}_multifrac.fixcol2.featureCounts",
    log:
        "logs/star/fct_fix_{aligner}_{trimmer}_{pese}_multifrac.log"
    threads: 1
    resources: time_min=220, mem_mb=8000, cpus=1
    shell: 
        "python3 scripts/fix_featurecounts_output.py -f {input} -c 2 > {output} 2> {log}"

