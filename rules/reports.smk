# AFTER MAPPING QUALITY CONTROL
#
# PER SAMPLE DNA
#
# https://multiqc.info/modules/bismark/

def qc_picard_DNA_input(wildcards):
    input = {}
    input["bam"] = "mapped/{sample}.bam"
    input["ref"] = expand("{ref_dir}/seq/{ref}.fa",ref_dir=reference_directory,ref=config["reference"])[0]
    if "lib_ROI" in config and config["lib_ROI"] != "wgs":
        input['lib_ROI'] = expand("{ref_dir}/intervals/{lib_ROI}/{lib_ROI}.interval_list",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0]
    return input

rule qc_picard_DNA:
    input:  unpack(qc_picard_DNA_input)
    output: table = "qc_reports/{sample}/qc_picard_DNA/picard.tsv",
    log:    "logs/{sample}/qc_picard_DNA.log"
    params: per_target = "qc_reports/{sample}/qc_picard_DNA/picard.per_target.tsv",
            wgs_chart = "qc_reports/{sample}/qc_picard_DNA/picard.wgs_chart.pdf",
            lib_ROI = config["lib_ROI"]
    threads: 1
    resources:  mem = 20
    conda: "../wrappers/qc_picard_DNA/env.yaml"
    script: "../wrappers/qc_picard_DNA/script.py"

def qc_qualimap_DNA_input(wildcards):
    input = {}
    input["bam"] = "mapped/{sample}.bam"
    if "lib_ROI" in config and config["lib_ROI"] != "wgs":
        input['lib_ROI'] = expand("{ref_dir}/intervals/{lib_ROI}/{lib_ROI}.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0]
    return input

rule qc_qualimap_DNA:
    input:  unpack(qc_qualimap_DNA_input)
    output: html = "qc_reports/{sample}/qc_qualimap_DNA/{sample}/qualimapReport.html"
    log:    "logs/{sample}/qc_qualimap_DNA.log"
    params: lib_ROI = config["lib_ROI"]
    threads: 4
    resources:  mem = 16
    conda: "../wrappers/qc_qualimap_DNA/env.yaml"
    script: "../wrappers/qc_qualimap_DNA/script.py"

rule qc_samtools:
    input:  bam = "mapped/{sample}.bam"
    output: idxstats = "qc_reports/{sample}/qc_samtools/{sample}.idxstats.tsv",
            flagstats = "qc_reports/{sample}/qc_samtools/{sample}.flagstat.tsv"
    log:    "logs/{sample}/qc_samtools.log"
    threads: 1
    conda: "../wrappers/qc_samtools/env.yaml"
    script: "../wrappers/qc_samtools/script.py"


###

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





# add input function with variable inputs... or make completely dependent on config["methylation_calling"]
# rule bismark2report_pe:
#     input:
#         alignment_report="results/reports/{sample}/{sample}_PE_report.txt",
#         nucleotide_report="results/reports/{sample}/{sample}.nucleotide_stats.txt",
#         mbias_report="methylation_calling/{sample}/{sample}.M-bias.txt",
#         splitting_report="methylation_calling/{sample}/{sample}_splitting_report.txt"
#     output:
#         html="methylation_calling/{sample}/{sample}.bismark2report.html"
#     log:
#         "logs/bismark/{sample}_bismark2report_no_deduplication.log"
#     params:
#         skip_optional_reports=True
#     conda:
#         "../wrappers/bismark/env.yaml"
#     script:
#         "../wrappers/bismark/bismark2report_pe.py"
