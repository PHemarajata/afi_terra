FROM mambaorg/micromamba:1.5.10-jammy

# fastp  → staphb/fastp    (FastpClean task)
# minimap2 + samtools are both included here:
#   - minimap2 for MinimapRick16S alignment
#   - samtools for BAM sort/index in MinimapRick16S and for
#     ExtractMetrics (Python subprocess calls to samtools idxstats/depth)
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    python=3.11 \
    pandas \
    minimap2 \
    samtools \
    && micromamba clean --all --yes

WORKDIR /opt/afi
COPY scripts /opt/afi/scripts

ENV PATH="/opt/conda/bin:/opt/afi/scripts:${PATH}"
