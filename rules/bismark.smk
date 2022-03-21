########################## BISMARK RULES

#copy reference file for BISMARK, cause BISMARK is stupid and cant put indexes anywhere else
rule copy_reference:
    input:
        expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
    output:
        expand("bismark/references/{ref_name}.fa",ref_name=config["reference"])[0],
    shell:
        "cp {input} {output}"

######################################################################################
# This script needs to be run only once to prepare the genome of interest for bisulfite alignments.
#     You need to specify a directory containing the genome you want to align your reads against
#     (please be aware that the bismark_genome_preparation script expects FastA files in this folder
#     (with either .fa or .fasta extension, single or multiple sequence entries per file).
#     Bismark will create two individual folders within this directory, one for a C->T converted genome and
#     the other one for the G->A converted genome. After creating C->T and G->A versions of the genome
#     they will be indexed in parallel using the indexer bowtie2-build (or hisat2-build). Once both C->T
#     and G->A genome indices have been created you do not need to use the genome preparation script again
#     (unless you want to align against a different genome...).
#
# Again, please note that Bowtie 2 and HISAT2 indexes are not compatible! To create a genome index for
#     use with HISAT2 the option --hisat2 needs to be included as well.

# likely will need - sudo apt-get install libtbb2

rule bismark_genome_preparation:
    input:
        expand("bismark/references/{ref_name}.fa",ref_name=config["reference"])[0],
    output:
        directory("bismark/references/Bisulfite_Genome")
    log:
        "logs/bismark/bismark_genome_preparation.log"
    params:
        "--verbose --bowtie2"  # optional params string
    threads: 
        6
    conda:  "../wrappers/bismark/env.yaml"
    script: "../wrappers/bismark/bismark_genome_preparation.py"

# The script bam2nuc reads BAM files and calculates the mono- and di-nucleotide coverage of the reads
# (using the genomic sequence rather than the observed sequence in the reads themselves) and
# compares it to the average genomic sequence composition. Reads harbouring InDels are not taken into consideration.
# Mono- or dinucleotides containing Ns are ignored as well.
rule bam2nuc_for_genome:
    input:
        genome_fa=expand("bismark/references/{ref_name}.fa",ref_name=config["reference"])[0],
    output:
        report="reports/genomic_nucleotide_frequencies.txt"
    log:
        "logs/bismark/bismark2nuc.log"
    threads: 
        6
    conda: "../wrappers/bismark/env.yaml"
    script: "../wrappers/bismark/bam2nuc_for_genome.py"

######################################################################################
# ALIGNMENT
# https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html#ii-bismark-alignment-step

def bismark_reads_input(wildcards):
    if config["is_paired"] == True:
        return {'fq_1': expand("trimmed_reads/{SAMPLE}_R1_val_1.fq.gz",SAMPLE=sample_tab.loc[sample_tab.sample_name == wildcards.sample_name, "sample_name"]),
                'fq_2': expand("trimmed_reads/{SAMPLE}_R2_val_1.fq.gz",SAMPLE=sample_tab.loc[sample_tab.sample_name == wildcards.sample_name, "sample_name"])}
    else:
        return {'fq': expand("trimmed_reads/{SAMPLE}_val_1.fq.gz",SAMPLE=sample_tab.loc[sample_tab.sample_name == wildcards.sample_name, "sample_name"])}

def bismark_alignment_report(wildcards):
    if config["is_paired"] == True:
        return expand("reports/{SAMPLE}/bismark/{SAMPLE}_PE_report.txt",SAMPLE=sample_tab.loc[sample_tab.sample_name == wildcards.sample_name, "sample_name"])
    else:
        return expand("reports/{SAMPLE}/bismark/{SAMPLE}_SE_report.txt",SAMPLE=sample_tab.loc[sample_tab.sample_name == wildcards.sample_name, "sample_name"])

rule bismark_alignment:
    input:
        reads=unpack(bismark_reads_input),
        genome=expand("bismark/references/{ref_name}.fa",ref_name=config["reference"])[0],
        bismark_indexes_dir="bismark/references/Bisulfite_Genome",
        genomic_freq="reports/genomic_nucleotide_frequencies.txt"
    output:
        bam=temp("bismark/{SAMPLE}/{SAMPLE}.bam"),
        report= bismark_alignment_report,
        nucleotide_stats="reports/{SAMPLE}/bismark/{SAMPLE}.nucleotide_stats.txt"
        # bam_unmapped_1="bismark/{SAMPLE}/{SAMPLE}_unmapped_reads_1.fq.gz",
        # bam_unmapped_2="bismark/{SAMPLE}/{SAMPLE}_unmapped_reads_2.fq.gz",
        # ambiguous_1="bismark/{SAMPLE}/{SAMPLE}_ambiguous_reads_1.fq.gz",
        # ambiguous_2="bismark/{SAMPLE}/{SAMPLE}_ambiguous_reads_2.fq.gz"
    log:
        "logs/bismark/{SAMPLE}_bismark_aligning.log"
    threads: 
        3
    params:
        extra="--nucleotide_coverage",
        basename="{SAMPLE}"
    conda: "../wrappers/bismark/env.yaml"
    script: "../wrappers/bismark/bismark_alignment.py"

######################################################################################
# DEDUPLICATION
# This command will deduplicate the Bismark alignment BAM file and remove all reads but one which align to the the very
# same position and in the same orientation. This step is recommended for whole-genome bisulfite samples,
# but should not be used for reduced representation libraries such as RRBS, amplicon or target enrichment libraries.

rule deduplicate_bismark:
    input: "bismark/{SAMPLE}/{SAMPLE}.bam"
    output:
        bam=temp("bismark/{SAMPLE}/{SAMPLE}.deduplicated.bam"),
        report="reports/{SAMPLE}/bismark/{SAMPLE}_deduplication_report.txt",
    log:
        "logs/bismark/{SAMPLE}_deduplication.log"
    params:
        extra=""  # optional params string
    threads: 
        6
    conda: "../wrappers/bismark/env.yaml"
    script: "../wrappers/bismark/bismark_deduplication.py"

rule bismark_methylation_extractor_dedup:
    input:
        "bismark/{SAMPLE}/{SAMPLE}.deduplicated.bam"
    output:
        #reports
        mbias_r1="reports/{SAMPLE}/bismark/{SAMPLE}.M-bias_R1.png",
        mbias_r2="reports/{SAMPLE}/bismark/{SAMPLE}.deduplicated.M-bias_R2.png",
        mbias_report="reports/{SAMPLE}/bismark/{SAMPLE}.deduplicated.M-bias.txt",
        splitting_report="reports/{SAMPLE}/bismark/{SAMPLE}.deduplicated.M-bias.txt",
        #methylation calls
        # coverafe file - 1-based chromosome coordinates ('inclusive') methylation info: % and counts
        methylome_CpG_cov="bismark/methylation_deduplicated/{SAMPLE}/{SAMPLE}.bismark.cov.gz",
            # BedGraph with methylation percentage: 0-based start, 1-based end (exclusive)
        methylome_CpG_mlevel_bedGraph="bismark/methylation_deduplicated/{SAMPLE}/{SAMPLE}.bedGraph.gz",
            # Primary output files: methylation status at each read cytosine position: (extremely large)
        read_base_meth_state_cpg="bismark/methylation_deduplicated/{SAMPLE}/CpG_context_{SAMPLE}.txt.gz",
            # * You could merge CHG, CHH using: --merge_non_CpG
        read_base_meth_state_chg="bismark/methylation_deduplicated/{SAMPLE}/CHG_context_{SAMPLE}.txt.gz",
        read_base_meth_state_chh="bismark/methylation_deduplicated/{SAMPLE}/CHH_context_{SAMPLE}.txt.gz"
        # genome wide cytosine report; contains strand information, ideal for bsses read.bismark(); 1-based chromosome coordinates 
        #cytosine_report="results/methylation/{SAMPLE}/{SAMPLE}.cytosine_report.txt.gz",
    log:
        "logs/bismark/{SAMPLE}_methylation_calling_deduplicated.log"
    threads: 
        2
    params:
        extra="--gzip --comprehensive --bedGraph --multicore 2"  # optional params string
    conda: "../wrappers/bismark/env.yaml"
    script: "../wrappers/bismark/bismark_methylation_extractor.py"

rule bismark2report_pe_dedup:
    input:    
        alignment_report=bismark_alignment_report,
        nucleotide_report="reports/{SAMPLE}/bismark/{SAMPLE}.nucleotide_stats.txt",
        mbias_report="reports/{SAMPLE}/bismark/{SAMPLE}.deduplicated.M-bias.txt",
        splitting_report="reports/{SAMPLE}/bismark/{SAMPLE}.deduplicated_splitting_report.txt",
    output:
        html="reports/{SAMPLE}/bismark/{SAMPLE}.deduplicated.bismark2report.html"
    log:
       "logs/bismark/{SAMPLE}_bismark2report_deduplicated.log"
    params:
        skip_optional_reports=True
    conda: "../wrappers/bismark/env.yaml"
    script: "../wrappers/bismark/bismark2report.py"

######################################################################################
# NO DEDUPLICATION - DEFAULT

rule bismark_methylation_extractor:
    input: 
        "bismark/{SAMPLE}/{SAMPLE}.bam"
    output:
        #reports
        mbias_r1="reports/{SAMPLE}/bismark/{SAMPLE}.M-bias_R1.png",
        mbias_r2="reports/{SAMPLE}/bismark/{SAMPLE}.M-bias_R2.png",
        mbias_report="reports/{SAMPLE}/bismark/{SAMPLE}.M-bias.txt",
        splitting_report="reports/{SAMPLE}/bismark/{SAMPLE}_splitting_report.txt",
        #methylation calls
            # coverafe file - 1-based chromosome coordinates ('inclusive') methylation info: % and counts        
        methylome_CpG_cov="bismark/methylation_calling/{SAMPLE}/{SAMPLE}.bismark.cov.gz",
            # BedGraph with methylation percentage: 0-based start, 1-based end (exclusive)
        methylome_CpG_mlevel_bedGraph="bismark/methylation_calling/{SAMPLE}/{SAMPLE}.bedGraph.gz",
            # Primary output files: methylation status at each read cytosine position: (extremely large)
        read_base_meth_state_cpg="bismark/methylation_calling/{SAMPLE}/CpG_context_{SAMPLE}.txt.gz",
            # * You could merge CHG, CHH using: --merge_non_CpG
        read_base_meth_state_chg="bismark/methylation_calling/{SAMPLE}/CHG_context_{SAMPLE}.txt.gz",
        read_base_meth_state_chh="bismark/methylation_calling/{SAMPLE}/CHH_context_{SAMPLE}.txt.gz"
        # genome wide cytosine report; contains strand information, ideal for bsses read.bismark(); 1-based chromosome coordinates 
        #cytosine_report="results/methylation/{SAMPLE}/{SAMPLE}.cytosine_report.txt.gz",
    log:
        "logs/bismark/{SAMPLE}_methylation_calling_no_deduplication.log"
    threads: 
        2
    params:
        extra="--gzip --comprehensive --bedGraph --multicore 2"  # optional params string
    conda: "../wrappers/bismark/env.yaml"
    script: "../wrappers/bismark/bismark_methylation_extractor.py"

rule bismark2report:
    input:    
        alignment_report=bismark_alignment_report,
        nucleotide_report="reports/{SAMPLE}/bismark/{SAMPLE}.nucleotide_stats.txt",
        mbias_report="reports/{SAMPLE}/bismark/{SAMPLE}.M-bias.txt",
        splitting_report="reports/{SAMPLE}/bismark/{SAMPLE}_splitting_report.txt",
    output:
        html="reports/{SAMPLE}/bismark/{SAMPLE}.bismark2report.html"
    log:
       "logs/bismark/{SAMPLE}_bismark2report_no_deduplication.log"
    params:
        skip_optional_reports=True
    conda: "../wrappers/bismark/env.yaml"
    script: "../wrappers/bismark/bismark2report.py"