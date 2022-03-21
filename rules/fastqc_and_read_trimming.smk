########################################################################################################################
# SINGLE END
rule fastqc_before_trimming_SE:
    input:
        "results/raw_reads/{SAMPLE}.fastq.gz",
    output:
        "results/reports/{SAMPLE}/fastqc_raw/{SAMPLE}_fastqc.html",
    log:
        "logs/{SAMPLE}/fastqc_before_trimming.log"
    threads:
        4
    params:
        FASTQCDIR="results/reports/{SAMPLE}/fastqc_raw"
    conda: "../envs/fastqc.yaml"
    shell:
        "fastqc --threads {threads} --quiet {input} -o {params.FASTQCDIR} >> {log} 2>&1"

# TRIMMING RAW READS - ADAPTERS + QUALITY <20
rule trim_galore_SE:
    input:
        "results/raw_reads/{SAMPLE}.fastq.gz",
    output:
        "results/trimmed_reads/{SAMPLE}_val_1.fq.gz",
        "results/trimmed_reads/{SAMPLE}.fastq.gz_trimming_report.txt",
    threads:
        2
    params:
        extra="--illumina --quality 20 --cores 2",
        rrbs=config["rrbs"],
        directional=config["rrbs_directional"],
    log:
        "logs/{SAMPLE}/trim_galore.log"
    conda: "../wrappers/trim_galore/env.yaml"
    script: "../wrappers/trim_galore/trim_galore_se.py"

# FASTQC TRIMMED READS
rule fastqc_after_trimming_SE:
    input:
        "results/trimmed_reads/{SAMPLE}_val_1.fq.gz",
    output:
        "results/reports/{SAMPLE}/fastqc_trimmed/{SAMPLE}_val_1_fastqc.html",
    log:
        "logs/{SAMPLE}/fastqc_after_trimming.log"
    threads:
        4
    params:
        FASTQCDIR="results/reports/{SAMPLE}/fastqc_trimmed"
    conda: "../envs/fastqc.yaml"
    shell:
        "fastqc --threads {threads} --quiet {input} -o {params.FASTQCDIR} >> {log} 2>&1"


########################################################################################################################
# PAIRED END
rule fastqc_before_trimming_PE:
    input:
        "results/raw_reads/{SAMPLE}_R1.fastq.gz",
        "results/raw_reads/{SAMPLE}_R2.fastq.gz"
    output:
        "results/reports/{SAMPLE}/fastqc_raw/{SAMPLE}_R1_fastqc.html",
        "results/reports/{SAMPLE}/fastqc_raw/{SAMPLE}_R2_fastqc.html"
    log:
        "logs/{SAMPLE}/fastqc_before_trimming.log"
    threads:
        4
    params:
        FASTQCDIR="results/reports/{SAMPLE}/fastqc_raw"
    conda: "../envs/fastqc.yaml"
    shell:
        "fastqc --threads {threads} --quiet {input} -o {params.FASTQCDIR} >> {log} 2>&1"

# TRIMMING RAW READS - ADAPTERS + QUALITY <20
rule trim_galore_PE:
    input:
        "results/raw_reads/{SAMPLE}_R1.fastq.gz",
        "results/raw_reads/{SAMPLE}_R2.fastq.gz"
    output:
        "results/trimmed_reads/{SAMPLE}_R1_val_1.fq.gz",
        "results/trimmed_reads/{SAMPLE}_R1.fastq.gz_trimming_report.txt",
        "results/trimmed_reads/{SAMPLE}_R2_val_2.fq.gz",
        "results/trimmed_reads/{SAMPLE}_R2.fastq.gz_trimming_report.txt"
    threads:
        2
    params:
        extra="--illumina --quality 20 --cores 2"
    log:
        "logs/{SAMPLE}/trim_galore.log"
    conda: "../wrappers/trim_galore/env.yaml"
    script: "../wrappers/trim_galore/trim_galore_pe.py"

# FASTQC TRIMMED READS
rule fastqc_after_trimming_PE:
    input:
        "results/trimmed_reads/{SAMPLE}_R1_val_1.fq.gz",
        "results/trimmed_reads/{SAMPLE}_R2_val_2.fq.gz"
    output:
        "results/reports/{SAMPLE}/fastqc_trimmed/{SAMPLE}_R1_val_1_fastqc.html",
        "results/reports/{SAMPLE}/fastqc_trimmed/{SAMPLE}_R2_val_2_fastqc.html"
    log:
        "logs/{SAMPLE}/fastqc_after_trimming.log"
    threads:
        4
    params:
        FASTQCDIR="results/reports/{SAMPLE}/fastqc_trimmed"
    conda: "../envs/fastqc.yaml"
    shell:
        "fastqc --threads {threads} --quiet {input} -o {params.FASTQCDIR} >> {log} 2>&1"

