# image: leviosam_env
FROM ubuntu:20.04
LABEL org.opencontainers.image.authors="naechyun@jhu.edu"

ENV TZ=America/New_York
ENV VERSION 0.4.2
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN apt-get update && apt-get install -y curl git build-essential cmake libsdsl-dev wget libbz2-dev zlib1g-dev liblzma-dev libncurses-dev zip
RUN wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2 && \
    tar -vxjf htslib-1.16.tar.bz2 && \
    cd htslib-1.16 && \
    ./configure && \
    make && \
    make install && \
    cd ../ && \
    rm htslib-1.16.tar.bz2
RUN wget https://github.com/samtools/samtools/releases/download/1.16/samtools-1.16.tar.bz2 && \
    tar -vxjf samtools-1.16.tar.bz2 && \
    cd samtools-1.16 && \
    ./configure && \
    make && \
    make install && \
    cd ../ && \
    rm samtools-1.16.tar.bz2
RUN wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 && \
    tar -vxjf bwa-0.7.17.tar.bz2 && \
    cd bwa-0.7.17 && \
    make && \
    cp bwa /usr/local/bin && \
    cd ../ && \
    rm bwa-0.7.17.tar.bz2
RUN wget https://github.com/BenLangmead/bowtie2/archive/refs/tags/v2.4.5.tar.gz && \
    tar -vxzf v2.4.5.tar.gz && \
    cd bowtie2-2.4.5 && \
    make && \
    cp bowtie2 bowtie2-build bowtie2-inspect bowtie2-build-s bowtie2-build-l bowtie2-align-s bowtie2-align-l bowtie2-inspect-s bowtie2-inspect-l /usr/local/bin/ && \
    cd ../ && \
    rm v2.4.5.tar.gz
RUN wget https://github.com/lh3/minimap2/archive/refs/tags/v2.24.tar.gz && \
    tar -vxzf v2.24.tar.gz && \
    cd minimap2-2.24 && \
    make && \
    cp minimap2 /usr/local/bin && \
    cd ../ && \
    rm v2.24.tar.gz
RUN wget https://github.com/marbl/Winnowmap/archive/refs/tags/v2.03.tar.gz && \
    tar -vxzf v2.03.tar.gz && \
    cd Winnowmap-2.03 && \
    make && \
    cp bin/winnowmap bin/meryl /usr/local/bin && \
    cd ../ && \
    rm v2.03.tar.gz
RUN curl -k -L https://github.com/milkschen/leviosam2/archive/refs/tags/v${VERSION}.tar.gz -o leviosam2-v${VERSION}.tar.gz && \
    tar -xzf leviosam2-v${VERSION}.tar.gz && \
    mv leviosam2-${VERSION} leviosam2 && \
    cd leviosam2 && \
    mkdir build && \
    cd build && \
    cmake -D CMAKE_BUILD_TYPE="release" .. && \
    echo "#ifndef VERSION\n#define VERSION \"${VERSION}\"\n#endif" > ../src/version.hpp && \
    make && \
    make install
ENV LD_LIBRARY_PATH="/usr/local/lib/:${LD_LIBRARY_PATH}"

