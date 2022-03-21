#########################

rule samtools_sort:
    input:
        "bismark/{SAMPLE}/{SAMPLE}.bam"
    output:
       protected("bismark/{SAMPLE}/{SAMPLE}.sorted.bam")
    threads:
        12  # Samtools takes additional threads through its option -@
    conda: "../envs/samtools.yaml"
    log:
        "logs/samtools/{SAMPLE}_samtools_sort.log"
    shell:
        """
        samtools sort -@ {threads} -o {output} -T {wildcards.SAMPLE} {input} >> {log} 2>&1
        # rm {input}
        """

rule samtools_index:
    input:
        "bismark/{SAMPLE}/{SAMPLE}.sorted.bam"
    output:
        protected("bismark/{SAMPLE}/{SAMPLE}.sorted.bam.bai")
    conda: "../envs/samtools.yaml"
    log:
        "logs/samtools/{SAMPLE}_samtools_index.log"
    shell:
       "samtools index {input} {output} >> {log} 2>&1"


rule samtools_sort_dedup:
    input:
        "bismark/{SAMPLE}/{SAMPLE}.deduplicated.bam"
    output:
        protected("bismark/{SAMPLE}/{SAMPLE}.deduplicated.sorted.bam")
    threads: 
        12 # Samtools takes additional threads through its option -@
    params:
        prefix="{SAMPLE}.deduplicated"
    conda: "../envs/samtools.yaml"
    log:
        "logs/samtools/{SAMPLE}_samtools_sort_dedup.log"
    shell:
        """
        samtools sort -@ {threads} -o {output} -T {params.prefix} {input} >> {log} 2>&1
        # rm {input}
        """

rule samtools_index_dedup:
    input:
        "bismark/{SAMPLE}/{SAMPLE}.deduplicated.sorted.bam"
    output:
        protected("bismark/{SAMPLE}/{SAMPLE}.deduplicated.sorted.bam.bai")
    conda: "../envs/samtools.yaml"
    log:
        "logs/samtools/{SAMPLE}_samtools_index_dedup.log"
    shell: 
       "samtools index {input} {output} >> {log} 2>&1"
