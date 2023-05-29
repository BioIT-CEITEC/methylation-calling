######################################
# wrapper for rule: bismark alignment
######################################
import re
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: bismakr alignment \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()





bwa_ref_prefix = re.sub(".bwt$","",snakemake.input.ref)

command = "bwa mem -t "+str(snakemake.threads)+\
        " -R '@RG\\tID:"+snakemake.params.entity_name+"_"+snakemake.wildcards.sample+"\\tSM:"+snakemake.wildcards.sample+"\\tPL:illumina' -v 1 " +\
        bwa_ref_prefix + " " + " ".join(snakemake.input.fastqs) + " 2>> " + log_filename + " | " +\
        "samtools sort -@ " +str(snakemake.threads)+" -o "+snakemake.output.bam+" /dev/stdin 2>> "+log_filename
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)


#######

# definice
extra = snakemake.params.get("extra", "")
cmdline_args = ["bismark {extra} --bowtie2"]

outdir = os.path.dirname(snakemake.output.bam)
if outdir:
    cmdline_args.append("--output_dir {outdir}")

genome_indexes_dir = os.path.dirname(snakemake.input.bismark_indexes_dir)
cmdline_args.append("{genome_indexes_dir}")


# basename
if snakemake.params.get("basename", None):
    cmdline_args.append("--basename {snakemake.params.basename:q}")
    basename = snakemake.params.basename
else:
    basename = None

#!!!!
# reads input
single_end_mode = snakemake.input.get("fq", None)
if single_end_mode:
    # for SE data, you only have to specify read1 input by -i or --in1, and
    # specify read1 output by -o or --out1.
    cmdline_args.append("--se {snakemake.input.fq:q}")
    mode_prefix = "se"
    if basename is None:
        basename = basename_without_ext(snakemake.input.fq)
else:
    # for PE data, you should also specify read2 input by -I or --in2, and
    # specify read2 output by -O or --out2.
    cmdline_args.append("-1 {snakemake.input.fq_1:q} -2 {snakemake.input.fq_2:q}")
    mode_prefix = "pe"

    if basename is None:
        # default basename
        basename = basename_without_ext(snakemake.input.fq_1) + "_bismark_bt2"

# log
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
cmdline_args.append("{log}")

# run
shell(" ".join(cmdline_args))

# Move outputs into proper position.
expected_2_actual_paths = [
    (
        snakemake.output.bam,
        os.path.join(
            outdir, "{}{}.bam".format(basename, "" if single_end_mode else "_pe")
        ),
    ),
    (
        snakemake.output.report,
        os.path.join(
            outdir,
            "{}_{}_report.txt".format(basename, "SE" if single_end_mode else "PE"),
        ),
    ),
    (
        snakemake.output.get("nucleotide_stats", None),
        os.path.join(
            outdir,
            "{}{}.nucleotide_stats.txt".format(
                basename, "" if single_end_mode else "_pe"
            ),
        ),
    ),
]
log_append = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)
for (exp_path, actual_path) in expected_2_actual_paths:
    if exp_path and (exp_path != actual_path):
        shell("mv {actual_path:q} {exp_path:q} {log_append}")