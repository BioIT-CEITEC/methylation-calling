import os
import snakemake.io
from multiprocessing import cpu_count
import pandas as pd
import json
from snakemake.utils import min_version

####################################
min_version("5.18.0")
configfile: "config.json"

GLOBAL_REF_PATH = config["globalResources"]
GLOBAL_TMPD_PATH = config["globalTmpdPath"]

os.makedirs(GLOBAL_TMPD_PATH, exist_ok=True)

# Reference processing
#
if config["lib_ROI"] != "wgs":
    # setting reference from lib_ROI
    f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","lib_ROI.json"))
    lib_ROI_dict = json.load(f)
    f.close()
    config["reference"] = [ref_name for ref_name in lib_ROI_dict.keys() if isinstance(lib_ROI_dict[ref_name],dict) and config["lib_ROI"] in lib_ROI_dict[ref_name].keys()][0]


# setting organism from reference
f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","reference2.json"),)
reference_dict = json.load(f)
f.close()
config["species_name"] = [organism_name for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]
config["organism"] = config["species_name"].split(" (")[0].lower().replace(" ","_")
if len(config["species_name"].split(" (")) > 1:
    config["species"] = config["species_name"].split(" (")[1].replace(")","")



##### Config processing #####
# Folders
#
reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])

# Samples
#
sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")

if not config["is_paired"]:
    read_pair_tags = [""]
    SEPEtag="SE"
else:
    read_pair_tags = ["_R1","_R2"]
    SEPEtag = "PE"

# CONFIGS
# enrichment_method - RRBS_MspI RRBS_MseI target_enrichment whole_genome
# type_seq pbat_seq em_seq BS_seq

#deduplication not for RRBS
# It is important to note that deduplication is not recommended for RRBS, amplicon or other target enrichment-type libraries!
config["deduplication"]=False
if config["enrichment_method"]=="RRBS_MspI" or config["enrichment_method"]=="RRBS_MseI" or config["enrichment_method"]=="target_enrichment":
    config["deduplication"]=False

if config["lib_ROI"]=="wgs" or config["lib_ROI"]=="WGS":
    config["deduplication"] = True
else:
    config["deduplication"] = False

# reports
config["qc_qualimap_DNA"]=True
config["qc_samtools"]=True
config["qc_picard_DNA"]=True

#empty methylation_calling folder
os.mkdir("methylation_calling")

wildcard_constraints:
    sample = "|".join(sample_tab.sample_name),
    read_pair_tag = "(_R.)?"


##### Target rules #####
def get_ruleall_output(wildcards):
    if config["methylation_calling"]:
        output_list=expand("methylation_calling/{sample}/{sample}.bismark.cov.gz",sample=sample_tab.sample_name)
    return output_list

rule all:
    input:  expand("mapped/{sample}.bam",sample = sample_tab.sample_name),
            expand("mapped/{sample}.bam.bai", sample = sample_tab.sample_name),
            get_ruleall_output,
            "qc_reports/multiqc.html"


##### Modules #####

include: "rules/trimming.smk"
include: "rules/bismark.smk"
include: "rules/bismark_deduplication.smk"
include: "rules/bismark_methylation.smk"
include: "rules/reports.smk"
