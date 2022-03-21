######################################
# wrapper for rule: bismark_genome_preparation
######################################

__author__ = "Roman Chernyatchik"
__copyright__ = "Copyright (c) 2019 JetBrains"
__email__ = "roman.chernyatchik@jetbrains.com"
__license__ = "MIT"

from os import path
from snakemake.shell import shell

input_dir = path.dirname(snakemake.input[0])

params_extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell("bismark_genome_preparation {params_extra} {input_dir} {log}")