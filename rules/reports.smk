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


########################################################################################################################
def multiqc_report_input(wildcards):
    input = {}
    if wildcards.sample != "all_samples":
        if config["qc_qualimap_DNA"]:
            input['qc_qualimap_DNA'] = "qc_reports/{sample}/qc_qualimap_DNA/{sample}/qualimapReport.html"
        if config["qc_samtools"]:
            input['qc_samtools'] = "qc_reports/{sample}/qc_samtools/{sample}.idxstats.tsv"
        if config["qc_picard_DNA"]:
            input['qc_picard_DNA'] = "qc_reports/{sample}/qc_picard_DNA/picard.tsv"

        if config["methylation_calling"]:
            input['mbias_report'] = "qc_reports/{sample}/bismark/m_bias/M-bias.txt",
            input['splitting_report'] = "qc_reports/{sample}/bismark/meth_extract/sample_splitting_report.txt"

        input['bismark_report'] = os.path.join("qc_reports/{sample}/bismark/align/",SEPEtag,"_report.txt"),
        input['bismark_nucleotide_stats'] = "qc_reports/{sample}/bismark/bam2nuc/nucleotide_stats.txt"

    else:
        input['per_sample_reports'] = expand("qc_reports/{sample}/single_sample_alignment_report.html",sample=sample_tab.sample_name)

    return input


rule multiqc_report:
    input:  unpack(multiqc_report_input)
    output: html="qc_reports/{sample}/multiqc.html"
    log:    "logs/{sample}/multiqc.log"
    params: multiqc_config = workflow.basedir+"/wrappers/multiqc_report/multiqc_config.txt",
            multiqc_path = "qc_reports/{sample}/"
    conda: "../wrappers/multiqc_report/env.yaml"
    script: "../wrappers/multiqc_report/script.py"

########################################################################################################################
def per_sample_alignment_report_input(wildcards):
    input = {}
    input['multiqc'] = "qc_reports/{sample}/multiqc.html"

    if config["qc_qualimap_DNA"]:
        input['qc_qualimap_DNA'] = "qc_reports/{sample}/qc_qualimap_DNA/{sample}/qualimapReport.html"
    if config["qc_samtools"]:
        input['qc_samtools'] = "qc_reports/{sample}/qc_samtools/{sample}.idxstats.tsv"
    if config["qc_picard_DNA"]:
        input['qc_picard_DNA'] = "qc_reports/{sample}/qc_picard_DNA/picard.tsv"

    if config["methylation_calling"]:
        input['mbias_report'] = "qc_reports/{sample}/bismark/m_bias/M-bias.txt",
        input['splitting_report'] = "qc_reports/{sample}/bismark/meth_extract/sample_splitting_report.txt"

    input['bismark_report'] = os.path.join("qc_reports/{sample}/bismark/align/",SEPEtag,"_report.txt"),
    input['bismark_nucleotide_stats'] = "qc_reports/{sample}/bismark/bam2nuc/nucleotide_stats.txt"

    return input

rule per_sample_alignment_report:
    input:  unpack(per_sample_alignment_report_input)
    output: sample_report = "qc_reports/{sample}/single_sample_alignment_report.html",
    params: sample_name = "{sample}",
            config = "./config.json",
            paired = paired,
    conda: "../wrappers/per_sample_alignment_report/env.yaml"
    script: "../wrappers/per_sample_alignment_report/script.Rmd"


