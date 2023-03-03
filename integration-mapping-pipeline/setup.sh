#!/bin/bash

# Get the current dir.
if [ -n "$BASH_VERSION" ]; then
    DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
elif [ -n "$ZSH_VERSION" ]; then
    DIR=${0:a:h}  # https://unix.stackexchange.com/a/115431
else
	echo "Error: Unknown shell; cannot determine path to GeneGraphDB local repository"
fi

export INTEGRATION_MAPPING_REPO_DIR="${DIR}"
DOCKER_IMAGE="integration_mapping"

# docker aliases
alias integration_mapping_docker_build="docker build -t $DOCKER_IMAGE $INTEGRATION_MAPPING_REPO_DIR/env"

integration_mapping_docker_run_func() {
    INTEGRATION_MAPPING_REPO_DIR=${1}
    WORKDIR=$(realpath ${2})

    docker run -it --rm \
	    -v ${INTEGRATION_MAPPING_REPO_DIR}/snakemake:/integration-mapping-pipeline/snakemake \
	    -v ${WORKDIR}:/integration-mapping-pipeline/WORKDIR \
	    integration_mapping
}

integration_mapping_docker_snakemake_func() {
    INTEGRATION_MAPPING_REPO_DIR=${1}
    WORKDIR=$(realpath ${2})
    THREADS=${3}

    docker run \
	    -v ${INTEGRATION_MAPPING_REPO_DIR}/snakemake:/integration-mapping-pipeline/snakemake \
	    -v ${WORKDIR}:/integration-mapping-pipeline/WORKDIR \
	    integration_mapping snakemake -j ${THREADS} --keep-going
}

alias integration_mapping_docker_run="integration_mapping_docker_run_func $INTEGRATION_MAPPING_REPO_DIR"
alias integration_mapping_docker_snakemake="integration_mapping_docker_snakemake_func $INTEGRATION_MAPPING_REPO_DIR"

