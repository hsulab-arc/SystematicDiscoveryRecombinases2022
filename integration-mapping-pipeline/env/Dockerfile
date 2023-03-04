FROM continuumio/miniconda3

RUN apt-get update && \
    apt-get install -y \
    g++ \
    gcc \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    make \
    zlib1g-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# DEBUG
# RUN apt-get -y install vim less htop
# DEBUG

RUN conda update -n base -c defaults conda && \
    conda install -y -c conda-forge python=3.8 mamba &&  \
    mamba install -y \
    -c conda-forge \
    -c bioconda \
    -c r \
    -c defaults \
    bowtie2 \
    bwa \
    cd-hit \
    emboss \
    fastp \
    fastqc \
    htslib \
    r \
    r-essentials \
    samtools \
    snakemake \
    && mamba clean --all -y && conda clean --all -y


RUN mamba install -y \
    -c conda-forge \
    -c bioconda \
    bioconductor-biostrings \
    bioconductor-genomicranges \
    bioconductor-genomicranges \
    bioconductor-helloranges \
    biopython \
    click \
    editdistance \
    networkx \
    pandas \
    pybedtools \
    pysam \
    r-bedr \
    r-ggseqlogo \
    scipy \
    tqdm \
    umi_tools \
    && mamba clean --all -y && conda clean --all -y

RUN pip install git+https://github.com/pdhsulab/MGEfinder-integration-mapping-pipeline.git

WORKDIR /integration-mapping-pipeline/snakemake
