# image: leviosam_env
FROM biocontainers/bwa:v0.7.17_cv1 AS bwa

FROM ubuntu:20.04
MAINTAINER naechyun@jhu.edu

ENV TZ=America/New_York
ENV VERSION 0.1.0
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN apt-get update && apt-get install -y curl git build-essential cmake libsdsl-dev wget libbz2-dev zlib1g-dev liblzma-dev libncurses-dev
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
RUN curl -k -L https://github.com/milkschen/leviosam2/archive/refs/tags/v${VERSION}.tar.gz -o leviosam2-v${VERSION}.tar.gz && \
    tar -xzf leviosam2-v${VERSION}.tar.gz && \
    cd leviosam2-${VERSION} && \
    mkdir build && \
    cd build && \
    cmake .. && \
    echo "#ifndef VERSION\n#define VERSION \"${VERSION}\"\n#endif" > ../src/version.hpp && \
    make && \
    make install
COPY --from=bwa /opt/conda/bin/bwa /usr/local/bin
ENV LD_LIBRARY_PATH="/usr/local/lib/:${LD_LIBRARY_PATH}"

