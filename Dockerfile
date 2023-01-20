# Copyright (C) 2022,2023 Genome Research Ltd.
### Load images
FROM quay.io/biocontainers/minimap2:2.17--hed695b0_3 AS minimap2
FROM quay.io/biocontainers/sniffles:1.0.12--h8b12597_1 AS sniffles
FROM quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0 AS bedtools
FROM bhklab/samtools-1.9.0:latest AS samtools

### Copy from images to local
WORKDIR /usr/local/bin
COPY --from=minimap2 /usr/local/bin/minimap2 ./minimap2
COPY --from=sniffles /usr/local/bin/sniffles ./sniffles
COPY --from=bedtools /usr/local/bin/bedtools ./bedtools

### Install Java and Nextflow
RUN apt-get update -qq -y && apt-get install -y curl wget && rm -rf /var/lib/apt/lists/* && apt-get clean
RUN wget https://download.oracle.com/java/18/archive/jdk-18.0.2_linux-x64_bin.tar.gz
RUN tar xzf jdk-18.0.2_linux-x64_bin.tar.gz && rm jdk-18.0.2_linux-x64_bin.tar.gz
ENV JAVA_HOME=/usr/local/bin/jdk-18.0.2/
RUN update-alternatives --install /usr/bin/java java ${JAVA_HOME%*/}/bin/java 20000
RUN update-alternatives --install /usr/bin/javac javac ${JAVA_HOME%*/}/bin/javac 20000
RUN java -version
RUN curl -s https://get.nextflow.io | bash

### Test
RUN samtools --version
RUN minimap2 --version
RUN sniffles --help
RUN bedtools --version
RUN nextflow -version

### Set up for pipeline work
WORKDIR /nextflow
