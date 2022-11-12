# image: leviosam2
FROM leviosam2_env:latest AS lev

FROM ubuntu:20.04
MAINTAINER naechyun@jhu.edu

ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN apt-get update && apt-get install -y perl python3
COPY --from=lev /leviosam2 /lev
COPY --from=lev /usr/local/bin/leviosam2 /usr/local/bin/
COPY --from=lev /usr/local/bin/bgzip /usr/local/bin/
COPY --from=lev /usr/local/bin/bwa /usr/local/bin/
COPY --from=lev /usr/local/bin/bowtie2 /usr/local/bin/
COPY --from=lev /usr/local/bin/bowtie2-build /usr/local/bin/
COPY --from=lev /usr/local/bin/bowtie2-inspect /usr/local/bin/
COPY --from=lev /usr/local/bin/bowtie2-build-s /usr/local/bin/
COPY --from=lev /usr/local/bin/bowtie2-build-l /usr/local/bin/
COPY --from=lev /usr/local/bin/bowtie2-align-s /usr/local/bin/
COPY --from=lev /usr/local/bin/bowtie2-align-l /usr/local/bin/
COPY --from=lev /usr/local/bin/bowtie2-inspect-s /usr/local/bin/
COPY --from=lev /usr/local/bin/bowtie2-inspect-l /usr/local/bin/
COPY --from=lev /usr/local/bin/minimap2 /usr/local/bin/
COPY --from=lev /usr/lib/gcc/x86_64-linux-gnu/9/libgomp.a /usr/local/lib/
COPY --from=lev /usr/lib/gcc/x86_64-linux-gnu/9/libgomp.so /usr/local/lib/
COPY --from=lev /usr/lib/x86_64-linux-gnu/libgomp.so.1 /usr/local/lib/
COPY --from=lev /usr/lib/x86_64-linux-gnu/libgomp.so.1.0.0 /usr/local/lib/
COPY --from=lev /usr/local/bin/winnowmap /usr/local/bin/
COPY --from=lev /usr/local/bin/meryl /usr/local/bin/
COPY --from=lev /usr/local/bin/samtools /usr/local/bin/
COPY --from=lev /usr/local/lib/libhts.a /usr/local/lib
COPY --from=lev /usr/local/lib/libhts.so /usr/local/lib
COPY --from=lev /usr/local/lib/libhts.so.1.16 /usr/local/lib
COPY --from=lev /usr/local/lib/libhts.so.3 /usr/local/lib
COPY --from=lev /usr/lib/x86_64-linux-gnu/libsdsl.a /usr/local/lib
COPY --from=lev /usr/lib/x86_64-linux-gnu/libsdsl.so /usr/local/lib
COPY --from=lev /usr/lib/x86_64-linux-gnu/libsdsl.so.2.1.0 /usr/local/lib
COPY --from=lev /usr/lib/x86_64-linux-gnu/libsdsl.so.3 /usr/local/lib
ENV LD_LIBRARY_PATH="/usr/local/lib/:${LD_LIBRARY_PATH}"

