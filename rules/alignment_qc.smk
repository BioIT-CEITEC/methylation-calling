#QUALIMAP
rule qualimap:
    input:
        "bismark/{SAMPLE}/{SAMPLE}.sorted.bam"
    output:
        "reports/{SAMPLE}/qualimap/qualimapReport.html"
    params:
        OUTDIRECTORY="reports/{SAMPLE}/qualimap"
    threads:
        6
    log:
        "logs/qualimap/{SAMPLE}_qualimap.log"
    conda: "../envs/qualimap.yaml"
    shell:
        "qualimap bamqc -bam {input}-ip -nt {threads} -os -outdir {params.OUTDIRECTORY} --java-mem-size=8G >> {log} 2>&1"

########################################################################################################################
# SINGLE END
rule multiqc_se:
    input:
        "results/reports/{SAMPLE}/{SAMPLE}_fastqc.html",
        "results/reports/{SAMPLE}/{SAMPLE}_val_1_fastqc.html",
        "results/reports/qualimap/{SAMPLE}/qualimapReport.html"
    output:
        directory("results/reports/multiqc_se/{SAMPLE}")
    params:
        basename="{SAMPLE}",

    threads:
        6
    log:
        "logs/multiqc/{SAMPLE}_multiqc.log"
    conda: "../envs/multiqc.yaml"
    shell:
        """
        mkdir {output}
        multiqc --force -o {output} -n {params.basename} {params.INDIRECTORY} >> {log} 2>&1
        """

########################################################################################################################
# PAIRED END
rule multiqc_pe:
    input:
        "results/reports/{SAMPLE}/{SAMPLE}_R1_fastqc.html",
        "results/reports/{SAMPLE}/{SAMPLE}_R2_fastqc.html",
        "results/reports/{SAMPLE}/{SAMPLE}_R1_val_1_fastqc.html",
        "results/reports/{SAMPLE}/{SAMPLE}_R2_val_2_fastqc.html",
        "results/reports/qualimap/{SAMPLE}/qualimapReport.html"
    output:
        directory("results/reports/multiqc_pe/{SAMPLE}")
    params:
        basename="{SAMPLE}",
    threads:
        6
    log:
        "logs/multiqc/{SAMPLE}_multiqc_pe.log"
    conda:
        os.path.join(workflow.basedir, "envs/multiqc.yaml")
    shell:
        """
        mkdir {output}
        multiqc --force -o {output} -n {params.basename} >> {log} 2>&1
        """