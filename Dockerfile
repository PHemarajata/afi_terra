FROM mambaorg/micromamba:1.5.10-jammy

RUN micromamba install -y -n base -c conda-forge -c bioconda \
    python=3.11 \
    pandas \
    fastp \
    minimap2 \
    samtools \
    kraken2 \
    && micromamba clean --all --yes

WORKDIR /opt/afi
COPY scripts /opt/afi/scripts

ENV PATH="/opt/afi/scripts:${PATH}"
