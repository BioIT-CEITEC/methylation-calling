#copy reference file for BISMARK, cause BISMARK is stupid and cant put indexes anywhere else
rule copy_bedgraphs:
    input:
        "results/bismark_methylation_not_deduplicated/{SAMPLE}/{SAMPLE}.bedGraph.gz",
    output:
        "results/mnp_klasifikace/{SAMPLE}.bedGraph.gz"
    log:
        "logs/copy/{SAMPLE}_copy_bedgraphs.log"
    shell:
        "cp {input} {output} 2> {log} "

# copy files which are mounted into docker (the folder results/mnp_klasifikace )
rule copy_mnp_scripts:
    input:
        os.path.join(workflow.basedir,"scripts/mnp_klasifikace/AgilentMethylseq_to_mnp_BASH_V3.R"),
        os.path.join(workflow.basedir,"scripts/mnp_klasifikace/10k_probes.Rdata")       
    output:
        "results/mnp_klasifikace/AgilentMethylseq_to_mnp_BASH_V3.R",
        "results/mnp_klasifikace/10k_probes.Rdata"
    log:
        "logs/copy/copy_mnp_scripts.log"
    shell:
        """
        cp {input[0]} {output[0]} 2> {log}
        cp {input[1]} {output[1]} 2> {log}
        """

rule bedtools_intersect:
    input: 
        "results/mnp_klasifikace/{SAMPLE}.bedGraph.gz"
    output: 
        "results/mnp_klasifikace/{SAMPLE}.intersected.bedGraph"
    params:
        probes10kfile=os.path.join(workflow.basedir,"scripts/mnp_klasifikace/10k_mnp_probes_extended_CpGislands.bed")
    conda: 
        os.path.join(workflow.basedir, "envs/bedtools.yaml")
    log:
        "logs/bedtools_intersect/{SAMPLE}.log"
    shell:
        """
        bedtools intersect -loj -a {params.probes10kfile} -b {input[0]} > {output[0]}  2> {log}
        """
# -loj Perform a “left outer join”. That is, for each feature in A report each overlap with B. If no overlaps are found, report a NULL feature for B.

rule mnp_klasifikace:
    input:
        "results/mnp_klasifikace/{SAMPLE}.bedGraph.gz",
        "results/mnp_klasifikace/{SAMPLE}.intersected.bedGraph",
        "results/mnp_klasifikace/AgilentMethylseq_to_mnp_BASH_V3.R",
        "results/mnp_klasifikace/10k_probes.Rdata"
    output:
        "results/mnp_klasifikace/{SAMPLE}_mnp_klasifikace.xlsx"
    params:
        mount = os.path.abspath("results/mnp_klasifikace/"),
        input = ["{SAMPLE}.bedGraph.gz","{SAMPLE}.intersected.bedGraph"]
    threads:
        12
    log:
        "logs/docker/{SAMPLE}_mnp_klasifikace_docker.log"
    script: 
        os.path.join(workflow.basedir, "scripts/mnp_klasifikace.py")
       