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
    threads: 10
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
    threads: 5
    conda: "../wrappers/qc_samtools/env.yaml"
    script: "../wrappers/qc_samtools/script.py"


########################################################################################################################
def multiqc_report_input(wildcards):
    input = {}
    input['qc_qualimap_DNA'] = expand("qc_reports/{sample}/qc_qualimap_DNA/{sample}/qualimapReport.html", sample = sample_tab.sample_name)
    input['qc_samtools'] = expand("qc_reports/{sample}/qc_samtools/{sample}.idxstats.tsv", sample = sample_tab.sample_name)
    input['qc_picard_DNA'] = expand("qc_reports/{sample}/qc_picard_DNA/picard.tsv", sample = sample_tab.sample_name)

    input['bismark_report'] = expand("qc_reports/{sample}/bismark/{sample}_{tag}_report.txt", sample = sample_tab.sample_name,tag=SEPEtag)
    input['bismark_nucleotide_stats'] = expand("qc_reports/{sample}/bismark/{sample}.nucleotide_stats.txt", sample = sample_tab.sample_name)

    if config["methylation_calling"]:
        input['mbias_report'] = expand("qc_reports/{sample}/bismark/{sample}_M-bias.txt", sample = sample_tab.sample_name)
        input['splitting_report'] = expand("qc_reports/{sample}/bismark/{sample}_splitting_report.txt", sample = sample_tab.sample_name)

    print(input)
    return input


rule multiqc_report:
    input:  unpack(multiqc_report_input)
    output: html="qc_reports/multiqc.html"
    log:    "logs/multiqc.log"
    params: multiqc_config = workflow.basedir+"/wrappers/multiqc_report/multiqc_config.txt",
        multiqc_path = "qc_reports/"
    conda: "../wrappers/multiqc_report/env.yaml"
    script: "../wrappers/multiqc_report/script.py"

########################################################################################################################


########################################################################################################################
#
# rule final_alignment_report:
#     input:  all_sample_multiqc = "qc_reports/multiqc.html"
#     output: html = "qc_reports/final_alignment_report.html"
#     conda: "../wrappers/final_alignment_report/env.yaml"
#     script: "../wrappers/final_alignment_report/script.Rmd"
