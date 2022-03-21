"""Snakemake wrapper for running docker container"""

__author__ = "Robin Jugas"
__copyright__ = "Copyright 2021, Robin Jugas"
__email__ = "robinjugas@gmail.com"
__license__ = "MIT"

import os
from snakemake.shell import shell



# Check that 3 input files were supplied
n = len(snakemake.input)
assert n == 4, "Input must contain 4 files. Given: %r." % n


# Check that four output files were supplied
m = len(snakemake.output)
assert m == 1, "Output must contain 1 file. Given: %r." % m

# Check that all input files are in the same directory
in_dir = os.path.dirname(snakemake.input[0])
for file_path in snakemake.input[1:]:
    assert in_dir == os.path.dirname(file_path), (
        "Input files need to be in one docker-mounted directory"
    )

shell(
    """
    if [ ! "$(docker ps -q -f name=MNP_DOCKER)" ]; then
        if [ "$(docker ps -aq -f status=exited -f name=MNP_DOCKER)" ]; then
            # cleanup
            docker stop MNP_DOCKER >> {snakemake.log} 2>&1
            docker rm MNP_DOCKER >> {snakemake.log} 2>&1
        fi
        # run your container
        docker run -d -e PASSWORD=1234 -v {snakemake.params.mount}:/home/data --name MNP_DOCKER mnp_klasifikace_docker:latest >> {snakemake.log} 2>&1
    fi
    """
)


shell(
    """
    docker exec --workdir /home/data MNP_DOCKER Rscript /home/data/AgilentMethylseq_to_mnp_BASH_V3.R {snakemake.params.input} >> {snakemake.log} 2>&1
    docker exec MNP_DOCKER chmod -R 777 /home/data >> {snakemake.log} 2>&1
    """
)

