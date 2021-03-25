# image: alshai/leviosam
FROM ubuntu:20.04
MAINTAINER tmun1@jhu.edu

ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN apt-get update && apt-get install -y git build-essential cmake libhts-dev libsdsl-dev
RUN git clone https://github.com/alshai/levioSAM.git /app
RUN mkdir -p /app/build && cd /app/build && cmake .. && make && make install
