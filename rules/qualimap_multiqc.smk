#QUALIMAP
rule qualimap:
    input:
        "results/bismark/{SAMPLE}/{SAMPLE}.sorted.bam"
    output:
        "results/reports/qualimap/{SAMPLE}/qualimapReport.html"
    params:
        BEDFILE=os.path.join(config["reference_folder"], config["bed_file"]),
        OUTDIRECTORY="results/reports/qualimap/{SAMPLE}"
    threads: 
        6
    log:
        "logs/qualimap/{SAMPLE}_qualimap.log"
    conda:
        os.path.join(workflow.basedir, "envs/qualimap.yaml")
    shell:
        "qualimap bamqc -bam {input} -gff {params.BEDFILE} -ip -nt {threads} -os -outdir {params.OUTDIRECTORY} --java-mem-size=4G >> {log} 2>&1"



rule multiqc:
    input:
        "results/reports/{SAMPLE}/{SAMPLE}_R1_fastqc.html",
        "results/reports/{SAMPLE}/{SAMPLE}_R2_fastqc.html",
        "results/reports/{SAMPLE}/{SAMPLE}_R1_val_1_fastqc.html",
        "results/reports/{SAMPLE}/{SAMPLE}_R2_val_2_fastqc.html",
        "results/reports/bismark_not_deduplicated/{SAMPLE}/{SAMPLE}.bismark2report.html",
        "results/reports/qualimap/{SAMPLE}/qualimapReport.html"
    output:
        directory("results/reports/multiqc/{SAMPLE}")
    params:
        basename="{SAMPLE}",
        INDIRECTORY="results/reports/bismark_not_deduplicated/{SAMPLE}/"
        # DIRECTORY1="results/reports/bismark_deduplicated/{SAMPLE}",        
    threads: 
        6
    log:
        "logs/multiqc/{SAMPLE}_multiqc.log"
    conda:
        os.path.join(workflow.basedir, "envs/multiqc.yaml")
    shell:
        """
        mkdir {output}
        multiqc --force -o {output} -n {params.basename} {params.INDIRECTORY} >> {log} 2>&1
        """


# rule qualimap_dedup:
#     input:
#         "results/bismark_deduplicated/{SAMPLE}/{SAMPLE}.deduplicated.sorted.bam"
#     output:
#         "results/reports/bismark_deduplicated/{SAMPLE}/qualimapReport.html"
#     params:
#         BEDFILE=os.path.join(config["reference_folder"], config["bed_file"]),
#         DIRECTORY="results/reports/bismark_deduplicated/{SAMPLE}/"
#     threads: 
#         6
#     log:
#         "logs/qualimap/{SAMPLE}_qualimap_dedup.log"
#     conda:
#         os.path.join(workflow.basedir, "envs/qualimap.yaml")
#     shell:
#         "qualimap bamqc -bam {input} -gff {params.BEDFILE} -ip -nt {threads} -os -outdir {params.DIRECTORY} --java-mem-size=4G >> {log} 2>&1"

