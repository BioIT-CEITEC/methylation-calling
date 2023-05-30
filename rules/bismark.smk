########################## REFERENCE

#copy reference file for BISMARK, cause BISMARK is stupid and cant put indexes anywhere else
rule copy_reference:
    input:
        expand("{ref_dir}/seq/{ref}.fa",ref_dir=reference_directory,ref=config["reference"])[0]
    output:
        expand("references/{ref}.fa",ref=config["reference"])[0]
    shell:
        "cp {input} {output}"

# likely will need - sudo apt-get install libtbb2
rule bismark_genome_preparation_fa_gz:
    input:
        expand("references/{ref}.fa",ref=config["reference"])[0]
    output:
        directory("references/Bisulfite_Genome")
    log:
        "logs/bismark/bismark_genome_preparation.log"
    params:
        "--verbose --bowtie2"  # optional params string
    threads: 
        10
    conda:
        "../wrappers/bismark/env.yaml"
    script: 
        "../wrappers/bismark/bismark_genome_preparation_fa_gz.py"

rule bam2nuc_for_genome:
    input:
        genome_fa=expand("references/{ref}.fa",ref=config["reference"])[0]
    output:
        report="references/genomic_nucleotide_frequencies.txt"
    log:
        "logs/bismark/bismark2nuc.log"
    threads: 
        10
    conda:
        "../wrappers/bismark/env.yaml"
    script: 
        "../wrappers/bismark/bam2nuc_for_genome.py"

########################## ALIGNMENT
def alignment_input(wildcards):
    if config["trim_adapters"]:
        preprocessed = "cleaned_fastq"
    else:
        preprocessed = "raw_fastq"

    if read_pair_tags == [""]:
        return {"fq":expand("{preprocessed}/{sample}.fastq.gz",preprocessed=preprocessed,sample=wildcards.sample)[0]}
    else:
        return {"fq_1":expand("{preprocessed}/{sample}_R1.fastq.gz",preprocessed=preprocessed,sample=wildcards.sample)[0],
        "fq_2":expand("{preprocessed}/{sample}_R2.fastq.gz",preprocessed=preprocessed,sample=wildcards.sample)[0]}


# CONFIGS
# enrichment_method - RRBS_MspI RRBS_MseI target_enrichment whole_genome
# type_seq pbat_seq em_seq BS_seq

bismark_alignment_params = [" --nucleotide_coverage -N 1 "]
if config["enrichment_method"] == "RRBS_MspI" or config["enrichment_method"] == "RRBS_MseI" :
    bismark_alignment_params.append(" --non_directional ")

if config["type_seq"]=="pbat_seq":
    bismark_alignment_params.append(" --pbat ")
elif config["type_seq"]=="em_seq":
    bismark_alignment_params.append(" --maxins 1000 ")


rule bismark_alignment:
    input:
        unpack(alignment_input),
        genome=expand("references/{ref}.fa",ref=config["reference"])[0],
        bismark_indexes_dir="references/Bisulfite_Genome",
        genomic_freq="references/genomic_nucleotide_frequencies.txt"
    output:
        bam="mapped/{sample}.not_markDups.bam",
        report="qc_reports/{sample}/bismark/{sample}_" + str(SEPEtag) + "_report.txt",
        nucleotide_stats="qc_reports/{sample}/bismark/{sample}.nucleotide_stats.txt"
    log:
        "logs/{sample}/bismark_alignment.log"
    threads: 
        5
    params:
        extra=bismark_alignment_params,
        basename="{sample}",
    conda:
        "../wrappers/bismark/env.yaml"
    script: 
        "../wrappers/bismark/bismark_alignment.py"







