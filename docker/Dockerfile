# image: alshai/leviosam
FROM ubuntu:20.04
MAINTAINER tmun1@jhu.edu

ENV TZ=America/New_York
ENV VERSION 0.4.0
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN apt-get update && apt-get install -y curl git build-essential cmake libhts-dev libsdsl-dev
RUN curl -k -L https://github.com/alshai/levioSAM/archive/refs/tags/v${VERSION}.tar.gz -o leviosam-v${VERSION}.tar.gz && \
    tar -xzf leviosam-v${VERSION}.tar.gz && \
    cd levioSAM-${VERSION} && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    make install
