import os
import snakemake.io
from multiprocessing import cpu_count


####################################

###########################################
# DEFINITION OF VARIABLES
# configfile: "config/config.yaml"

# SAMPLE=config["samples"]

R1R2 = ["_R1", "_R2"]

####################################
# RULES
include: "rules/fastqc_and_read_trimming.smk"
include: "rules/bismark.smk"
include: "rules/qualimap_multiqc.smk"
include: "rules/samtools.smk"


# RULE ALL
rule all:
    input:
        expand(["results/reports/{SAMPLE}/{SAMPLE}_R1_fastqc.html",
                "results/reports/{SAMPLE}/{SAMPLE}_R2_fastqc.html",
                "results/reports/{SAMPLE}/{SAMPLE}_R1_val_1_fastqc.html",
                "results/reports/{SAMPLE}/{SAMPLE}_R2_val_2_fastqc.html",

                # bismark not duplicated
                "results/bismark_methylation_not_deduplicated/{SAMPLE}/{SAMPLE}.bismark.cov.gz",
                "results/bismark_methylation_not_deduplicated/{SAMPLE}/{SAMPLE}.bedGraph.gz",
                "results/reports/bismark_not_deduplicated/{SAMPLE}/{SAMPLE}.bismark2report.html",

                # bismark deduplicated
                # "results/bismark_methylation_deduplicated/{SAMPLE}/{SAMPLE}.bismark.cov.gz",
                # "results/bismark_methylation_deduplicated/{SAMPLE}/{SAMPLE}.bedGraph.gz",                
                # "results/reports/bismark_deduplicated/{SAMPLE}/{SAMPLE}.bismark2report.html",
                # samtools sorting
                "results/bismark/{SAMPLE}/{SAMPLE}.sorted.bam",
                "results/bismark/{SAMPLE}/{SAMPLE}.sorted.bam.bai",
                # "results/bismark_deduplicated/{SAMPLE}/{SAMPLE}.deduplicated.sorted.bam",
                # "results/bismark_deduplicated/{SAMPLE}/{SAMPLE}.deduplicated.sorted.bam.bai",
                #qualimap  & reports
                "results/reports/qualimap/{SAMPLE}/qualimapReport.html",
                "results/reports/multiqc/{SAMPLE}/",

                ],
               SAMPLE=config["samples"])
