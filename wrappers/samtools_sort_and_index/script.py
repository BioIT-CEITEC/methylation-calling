######################################
# wrapper for rule: index_and_stats
######################################
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: index_and_stats \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()


command = "samtools sort -n -@ " + str(snakemake.threads) + " -o " +snakemake.output.bam + " " + str(snakemake.input)
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "samtools index -@ " + str(snakemake.threads) + " " + snakemake.output.bam
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

#
# command = "samtools idxstats "+snakemake.output.bam+" > "+snakemake.output.idxstats+" 2>> "+log_filename+" "
# f = open(log_filename, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)
#
# command ="samtools flagstat "+snakemake.output.bam+" > "+snakemake.output.flagstats+" 2>> "+log_filename+" "
# f = open(log_filename, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)


command = "rm " + str(snakemake.input)
f = open(log_filename, 'at')
f.write("## COMMAND: " + command + "\n")
f.close()
shell(command)