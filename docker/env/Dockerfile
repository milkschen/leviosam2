# image: leviosam_env
#FROM biocontainers/bwa:v0.7.17_cv1 AS bwa

FROM ubuntu:20.04
MAINTAINER naechyun@jhu.edu

ENV TZ=America/New_York
ENV VERSION 0.2.0
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN apt-get update && apt-get install -y curl git build-essential cmake libsdsl-dev wget libbz2-dev zlib1g-dev liblzma-dev libncurses-dev zip
RUN wget https://github.com/samtools/htslib/releases/download/1.13/htslib-1.13.tar.bz2 && \
    tar -vxjf htslib-1.13.tar.bz2 && \
    cd htslib-1.13 && \
    ./configure && \
    make && \
    make install && \
    cd ../
RUN wget https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2 && \
    tar -vxjf samtools-1.13.tar.bz2 && \
    cd samtools-1.13 && \
    ./configure && \
    make && \
    make install && \
    cd ../
RUN wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 && \
    tar -vxjf bwa-0.7.17.tar.bz2 && \
    cd bwa-0.7.17 && \
    make && \
    make install && cd \
    ../
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.4.5/bowtie2-2.4.5-linux-x86_64.zip && \
    unzip bowtie2-2.4.5-linux-x86_64.zip && \
    cd bowtie2-2.4.5-linux-x86_64 && \
    mv bowtie2 bowtie2-build bowtie2-inspect bowtie2-build-s bowtie2-build-l bowtie2-align-s bowtie2-align-l bowtie2-inspect-s bowtie2-inspect-l /usr/local/bin/ && \
    cd ../
RUN curl -k -L https://github.com/milkschen/leviosam2/archive/refs/tags/v${VERSION}.tar.gz -o leviosam2-v${VERSION}.tar.gz && \
    tar -xzf leviosam2-v${VERSION}.tar.gz && \
    cd leviosam2-${VERSION} && \
    mkdir build && \
    cd build && \
    cmake .. && \
    echo "#ifndef VERSION\n#define VERSION \"${VERSION}\"\n#endif" > ../src/version.hpp && \
    make && \
    make install
#COPY --from=bwa /opt/conda/bin/bwa /usr/local/bin
ENV LD_LIBRARY_PATH="/usr/local/lib/:${LD_LIBRARY_PATH}"

