# TRIMMING trim galore
def trim_adapters_input(wildcards):
    if config["trim_adapters"]:
        if read_pair_tags == [""]:
            return "raw_fastq/{sample}.fastq.gz"
        else:
            return ["raw_fastq/{sample}_R1.fastq.gz","raw_fastq/{sample}_R2.fastq.gz"]

# https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md
# Note for RRBS using MseI:
# If your DNA material was digested with MseI (recognition motif: TTAA)
# instead of MspI it is NOT necessary to specify --rrbs or --non_directional since
# virtually all reads should start with the sequence TAA, and this holds true for both directional
# and non-directional libraries. As the end-repair of TAA restricted sites does not involve
# any cytosines it does not need to be treated especially. Instead, simply run Trim Galore!
# in the standard, i.e. non-RRBS, mode.

# CONFIGS
# enrichment_method - RRBS_MspI RRBS_MseI target_enrichment whole_genome
# type_seq pbat_seq em_seq BS_seq


trim_galore_params = [""]
if config["enrichment_method"]=="RRBS_MspI":
    trim_galore_params.append(" --rrbs --non_directional ")

if config["type_seq"]=="pbat_seq":
    trim_galore_params.append(" --clip_R1 6 --clip_R2 9 --three_prime_clip_r1 6 --three_prime_clip_r2 9 ")
if config["type_seq"]=="em_seq":
    trim_galore_params.append(" -clip_R1 8 --clip_R2 8 --three_prime_clip_r1 8 --three_prime_clip_r2 8 ")


rule trim_adapters:
    input:  trim_adapters_input,
    output: fastq = expand("cleaned_fastq/{{sample}}{read_pair_tag}.fastq.gz", read_pair_tag = read_pair_tags),
            trim_stats = expand("qc_reports/{{sample}}/trim_galore/trim_stats{read_pair_tag}.log", read_pair_tag=read_pair_tags)
    log:    "logs/{sample}/trim_adapters.log"
    params: paired = config["is_paired"],
            outdir = "cleaned_fastq",
            parameters = trim_galore_params
    threads: 8
    conda: "../wrappers/trim_adapters/env.yaml"
    script: "../wrappers/trim_adapters/script.py"

