# image: alshai/leviosam
FROM ubuntu:20.04
MAINTAINER tmun1@jhu.edu

ENV TZ=America/New_York
ENV VERSION 0.5.1
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN apt-get update && apt-get install -y curl git build-essential cmake libsdsl-dev wget libbz2-dev zlib1g-dev liblzma-dev
RUN wget https://github.com/samtools/htslib/releases/download/1.13/htslib-1.13.tar.bz2 && \
    tar -vxjf htslib-1.13.tar.bz2 && \
    cd htslib-1.13 && \
    ./configure && \
    make && \
    make install && \
    cd ../
RUN curl -k -L https://github.com/alshai/levioSAM/archive/refs/tags/v${VERSION}.tar.gz -o leviosam-v${VERSION}.tar.gz && \
    tar -xzf leviosam-v${VERSION}.tar.gz && \
    cd levioSAM-${VERSION} && \
    mkdir build && \
    cd build && \
    cmake .. && \
    echo "#ifndef VERSION\n#define VERSION \"${VERSION}\"\n#endif" > ../src/version.hpp && \
    make && \
    make install
ENV LD_LIBRARY_PATH="/usr/local/lib/:${LD_LIBRARY_PATH}"
