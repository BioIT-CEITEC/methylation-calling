######################################
# wrapper for rule: trim_adapters
######################################
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: trim_adapters \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

if snakemake.params.paired:
        paired_flag = " --paired"
        r1 = snakemake.output[0]
        r2 = snakemake.output[1]
else:
        paired_flag = ""
        r1 = snakemake.output[0]

command = "trim_galore --fastqc " + paired_flag + " " + str(snakemake.input) +" -o "+snakemake.params.outdir+ " 2>> "+log_filename
with open(log_filename, 'at') as f:
        f.write("## COMMAND: " + command + "\n")
shell(command)


if snakemake.params.paired:
        command = "mv " + r1.replace(".fastq.gz","_val_1.fq.gz") + " " + r1 + " 2>> "+log_filename
        with open(log_filename, 'at') as f:
                f.write("## COMMAND: " + command + "\n")
        shell(command)

        command = "mv " + r2.replace(".fastq.gz","_val_2.fq.gz") + " " + r2 + " 2>> "+log_filename
        with open(log_filename, 'at') as f:
                f.write("## COMMAND: " + command + "\n")
        shell(command)

        # report for multiQC
        command = "cp " + r1+"_trimming_report.txt" + " " + snakemake.output.trim_stats[0] + " 2>> "+log_filename
        with open(log_filename, 'at') as f:
                f.write("## COMMAND: " + command + "\n")
        shell(command)

        command = "cp " + r2+"_trimming_report.txt" + " " + snakemake.output.trim_stats[1] + " 2>> "+log_filename
        with open(log_filename, 'at') as f:
                f.write("## COMMAND: " + command + "\n")
        shell(command)

else:
        command = "mv " + r1.replace(".fastq.gz","_trimmed.fq.gz") + " " + r1 + " 2>> "+log_filename
        with open(log_filename, 'at') as f:
                f.write("## COMMAND: " + command + "\n")
        shell(command)

        # report for multiQC
        command = "cp " + r1+"_trimming_report.txt" + " " + snakemake.output.trim_stats[0] + " 2>> "+log_filename
        with open(log_filename, 'at') as f:
                f.write("## COMMAND: " + command + "\n")
        shell(command)

