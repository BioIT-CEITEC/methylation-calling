#########################

rule samtools_sort:
    input:
        "results/bismark/{SAMPLE}/{SAMPLE}.bam"
    output:
       protected("results/bismark/{SAMPLE}/{SAMPLE}.sorted.bam")
    threads:
        12  # Samtools takes additional threads through its option -@
    conda:
        os.path.join(workflow.basedir, "envs/samtools.yaml")
    log:
        "logs/samtools/{SAMPLE}_samtools_sort.log"
    shell:
        """
        samtools sort -@ {threads} -o {output} -T {wildcards.SAMPLE} {input} >> {log} 2>&1
        # rm {input}
        """

rule samtools_index:
    input:
        "results/bismark/{SAMPLE}/{SAMPLE}.sorted.bam"
    output:
        protected("results/bismark/{SAMPLE}/{SAMPLE}.sorted.bam.bai")
    conda:
        os.path.join(workflow.basedir, "envs/samtools.yaml")
    log:
        "logs/samtools/{SAMPLE}_samtools_index.log"
    shell:
       "samtools index {input} {output} >> {log} 2>&1"


rule samtools_sort_dedup:
    input:
        "results/bismark_deduplicated/{SAMPLE}/{SAMPLE}.deduplicated.bam"
    output:
        protected("results/bismark_deduplicated/{SAMPLE}/{SAMPLE}.deduplicated.sorted.bam")
    threads: 
        12 # Samtools takes additional threads through its option -@
    params:
        prefix="{SAMPLE}.deduplicated"
    conda:
        os.path.join(workflow.basedir, "envs/samtools.yaml")
    log:
        "logs/samtools/{SAMPLE}_samtools_sort_dedup.log"
    shell:
        """
        samtools sort -@ {threads} -o {output} -T {params.prefix} {input} >> {log} 2>&1
        # rm {input}
        """

rule samtools_index_dedup:
    input:
        "results/bismark_deduplicated/{SAMPLE}/{SAMPLE}.deduplicated.sorted.bam"
    output:
        protected("results/bismark_deduplicated/{SAMPLE}/{SAMPLE}.deduplicated.sorted.bam.bai")
    conda:
        os.path.join(workflow.basedir, "envs/samtools.yaml")
    log:
        "logs/samtools/{SAMPLE}_samtools_index_dedup.log"
    shell: 
       "samtools index {input} {output} >> {log} 2>&1"
