# It is important to note that deduplication is not recommended for RRBS, amplicon or other target enrichment-type libraries!

#bismark deduplication rule
if config["deduplication"]:
    rule deduplicate_bismark:
        input: "mapped/{sample}.not_markDups.bam",
        output:
            bam="mapped/{sample}.markDups.bam",
            report="mapped/{sample}.deduplication_report.txt",
        log:
            "logs/{sample}/bismark_deduplication.log"
        params:
            extra=""  # optional params string
        threads:
            6
        conda:
            "../wrappers/bismark/env.yaml"
        script:
            "../wrappers/bismark/bismark_deduplication.py"

# samtools sort and index
def rename_input(wildcards):
    if config["deduplication"]:
        return "mapped/{sample}.markDups.bam"
    else:
        return "mapped/{sample}.not_markDups.bam"

rule rename:
    input:  rename_input
    output: "mapped/{sample}.bam"
    threads:    2
    shell:  "mv {input} {output}"

# rule samtools_sort_and_index: #not recommended cause methylation extractor requires unsorted BAM or sorted by reads name
#     input:  samtools_sort_and_index_input,
#     output: bam = "mapped/{sample}.bam",
#             bai = "mapped/{sample}.bam.bai",
#     log:    "logs/{sample}/samtools.log"
#     threads:    8
#     conda: "../wrappers/samtools_sort_and_index/env.yaml"
#     script: "../wrappers/samtools_sort_and_index/script.py"