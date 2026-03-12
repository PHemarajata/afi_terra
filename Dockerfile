FROM ubuntu:22.04

RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    samtools \
    minimap2 \
    fastp \
    wget \
    curl \
    git

RUN pip3 install pandas

WORKDIR /opt/afi
COPY scripts /opt/afi/scripts

ENV PATH="/opt/afi/scripts:$PATH"
