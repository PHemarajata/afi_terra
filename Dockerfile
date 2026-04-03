FROM mambaorg/micromamba:1.5.10-jammy

# fastp  → staphb/fastp    (FastpClean task)
# minimap2 → staphb/minimap2 (MinimapRick16S task, already bundles samtools)
# samtools is kept here for ExtractMetrics (Python subprocess calls to
#   samtools idxstats / samtools depth on the 16S alignment BAM).
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    python=3.11 \
    pandas \
    samtools \
    && micromamba clean --all --yes

WORKDIR /opt/afi
COPY scripts /opt/afi/scripts

ENV PATH="/opt/conda/bin:/opt/afi/scripts:${PATH}"
