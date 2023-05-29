# Methylation calling

# CONFIGS
# enrichment_method - RRBS_MspI RRBS_MseI target_enrichment whole_genome
# directionality - directional non_directional
# type_seq pbat_seq em_seq


bismark_methylation_params = [" --gzip --comprehensive --bedGraph --multicore 2 "]
if config["enrichment_method"]=="RRBS_MspI" or config["enrichment_method"]=="RRBS_MseI":
    bismark_methylation_params.append(" --ignore 5 --ignore_r2 5 --ignore_3prime 5 --ignore_3prime_r2 5 ")
# elif config["type_seq"]=="pbat_seq" or config["type_seq"]=="em_seq": #uz otrimovano v readech
#     bismark_methylation_params.append(" --ignore 5 --ignore_r2 5 --ignore_3prime 5 --ignore_3prime_r2 5 ")
else:
    bismark_methylation_params.append(" –-ignore_r2 2 --ignore_3prime_r2 2 ") #default


if read_pair_tags == [""]: #SINGLE END
    rule bismark_methylation_extractor_se:
        input:
            "mapped/{sample}.bam"
        output:
            mbias_r1="qc_reports/{sample}/bismark/m_bias/M-bias_R1.png",
            mbias_report="qc_reports/{sample}/bismark/m_bias/M-bias.txt",
            splitting_report="qc_reports/{sample}/bismark/meth_extract/sample_splitting_report.txt",
            # coverage file - 1-based chromosome coordinates ('inclusive') methylation info: % and counts
            methylome_CpG_cov="methylation_calling/{sample}/{sample}.bismark.cov.gz",
            # BedGraph with methylation percentage: 0-based start, 1-based end (exclusive)
            methylome_CpG_mlevel_bedGraph="methylation_calling/{sample}/{sample}.bedGraph.gz",
            # Primary output files: methylation status at each read cytosine position: (extremely large)
            read_base_meth_state_cpg="methylation_calling/{sample}/CpG_context_{sample}.txt.gz",
            # * You could merge CHG, CHH using: --merge_non_CpG
            read_base_meth_state_chg="methylation_calling/{sample}/CHG_context_{sample}.txt.gz",
            read_base_meth_state_chh="methylation_calling/{sample}/CHH_context_{sample}.txt.gz",
            # genome wide cytosine report; contains strand information, ideal for bsses read.bismark(); 1-based chromosome coordinates
            #cytosine_report="results/methylation/{sample}/{sample}.cytosine_report.txt.gz",
        log:
            "logs/{sample}/bismark_methylation_extraction_se.log"
        threads:
            10
        params:
            extra=bismark_methylation_params  # optional params string
        conda:
            "../wrappers/bismark/env.yaml"
        script:
            "../wrappers/bismark/bismark_methylation_extractor.py"

else:  # PAIRED END
    rule bismark_methylation_extractor_pe:
        input: "mapped/{sample}.bam"
    output:
        mbias_r1="qc_reports/{sample}/bismark/m_bias/M-bias_R1.png",
        mbias_r2="qc_reports/{sample}/bismark/m_bias/M-bias_R2.png",
        mbias_report="qc_reports/{sample}/bismark/m_bias/M-bias.txt",
        splitting_report="qc_reports/{sample}/bismark/meth_extract/sample_splitting_report.txt",
        # coverage file - 1-based chromosome coordinates ('inclusive') methylation info: % and counts
        methylome_CpG_cov="methylation_calling/{sample}/{sample}.bismark.cov.gz",
        # BedGraph with methylation percentage: 0-based start, 1-based end (exclusive)
        methylome_CpG_mlevel_bedGraph="methylation_calling/{sample}/{sample}.bedGraph.gz",
        # Primary output files: methylation status at each read cytosine position: (extremely large)
        read_base_meth_state_cpg="methylation_calling/{sample}/CpG_context_{sample}.txt.gz",
        # * You could merge CHG, CHH using: --merge_non_CpG
        read_base_meth_state_chg="methylation_calling/{sample}/CHG_context_{sample}.txt.gz",
        read_base_meth_state_chh="methylation_calling/{sample}/CHH_context_{sample}.txt.gz",
        # genome wide cytosine report; contains strand information, ideal for bsses read.bismark(); 1-based chromosome coordinates
        #cytosine_report="results/methylation/{sample}/{sample}.cytosine_report.txt.gz",
    log:
        "logs/{sample}/bismark_methylation_extraction_pe.log"
    threads:
        10
    params:
        extra=bismark_methylation_params  # optional params string
    conda:
        "../wrappers/bismark/env.yaml"
    script:
        "../wrappers/bismark/bismark_methylation_extractor.py"


# https://sequencing.qcfail.com/articles/library-end-repair-reaction-introduces-methylation-biases-in-paired-end-pe-bisulfite-seq-applications/
# Since this is a common phenomenon with any standard directional library the option –ignore_r2 2 of the methylation extractor should be added to the default data processing pipeline for paired-end applications.
# https: // sequencing.qcfail.com / articles / mispriming - in -pbat - libraries - causes - methylation - bias - and -poor - mapping - efficiencies /