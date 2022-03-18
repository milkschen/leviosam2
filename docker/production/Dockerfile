# image: leviosam2
FROM biocontainers/bwa:v0.7.17_cv1 AS bwa
FROM leviosam2_env:latest AS lev

FROM ubuntu:20.04
MAINTAINER naechyun@jhu.edu

ENV TZ=America/New_York
COPY --from=bwa /opt/conda/bin/bwa /usr/local/bin
COPY --from=lev /leviosam2-0.2.0 /lev
COPY --from=lev /usr/local/bin/leviosam2 /usr/local/bin/
COPY --from=lev /usr/local/bin/bgzip /usr/local/bin/
COPY --from=lev /usr/local/bin/bwa /usr/local/bin/
COPY --from=lev /usr/local/bin/samtools /usr/local/bin/
COPY --from=lev /usr/local/lib/libhts.a /usr/local/lib
COPY --from=lev /usr/local/lib/libhts.so /usr/local/lib
COPY --from=lev /usr/local/lib/libhts.so.1.13 /usr/local/lib
COPY --from=lev /usr/local/lib/libhts.so.3 /usr/local/lib
COPY --from=lev /usr/lib/x86_64-linux-gnu/libsdsl.a /usr/local/lib
COPY --from=lev /usr/lib/x86_64-linux-gnu/libsdsl.so /usr/local/lib
COPY --from=lev /usr/lib/x86_64-linux-gnu/libsdsl.so.2.1.0 /usr/local/lib
COPY --from=lev /usr/lib/x86_64-linux-gnu/libsdsl.so.3 /usr/local/lib
ENV LD_LIBRARY_PATH="/usr/local/lib/:${LD_LIBRARY_PATH}"

