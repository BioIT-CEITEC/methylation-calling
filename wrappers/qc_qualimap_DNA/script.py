######################################
# wrapper for rule: qc_qualimap_DNA
######################################
import os
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: qc_qualimap_DNA \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

command = "mkdir -p $(dirname "+snakemake.output.html+") >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "bash -c unset DISPLAY >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

if snakemake.params.lib_ROI != "wgs":
    lib_ROI_param = " -gff " +snakemake.input.lib_ROI
else:
    lib_ROI_param = ""

command = "qualimap bamqc -bam " + snakemake.input.bam + lib_ROI_param + " -outdir "+ os.path.dirname(snakemake.output.html) +" -nt "+str(snakemake.threads)+" \
    --java-mem-size="+str(snakemake.resources.mem)+"G >> "+log_filename+" 2>&1"

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
