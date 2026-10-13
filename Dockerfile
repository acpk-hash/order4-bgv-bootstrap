# Asiacrypt 2026 artifact, paper #528 (order-four BGV digit extraction).
# Toolchain matches the measurement host: Ubuntu 22.04, GCC 11.4,
# NTL 11.5.1 (built from source, threads + GMP), GMP 6.2.1.
#
#   docker build -t order4-bgv .
#   docker run --rm -it order4-bgv            # shell in /artifact
#   docker run --rm order4-bgv ae/run.sh smoke
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential cmake git ca-certificates wget m4 patch \
        libgmp-dev python3 python3-pip time procps \
    && rm -rf /var/lib/apt/lists/*

ARG NTL_VERSION=11.5.1
RUN wget -q https://libntl.org/ntl-${NTL_VERSION}.tar.gz \
    && tar xzf ntl-${NTL_VERSION}.tar.gz \
    && cd ntl-${NTL_VERSION}/src \
    && ./configure SHARED=on \
    && make -j"$(nproc)" && make install \
    && cd / && rm -rf ntl-${NTL_VERSION}* && ldconfig
# NTL defaults, as on the measurement host. In particular NTL_ENABLE_AVX_FFT
# must stay off: it shrinks NTL's single-precision primes from 60 to 50 bits,
# which changes the modulus chain HElib builds (1066 -> 1086 bits on m=2^16).
COPY ae/check_ntl.cpp /tmp/check_ntl.cpp
RUN g++ /tmp/check_ntl.cpp -o /tmp/check_ntl -lntl -lgmp -lpthread && /tmp/check_ntl

WORKDIR /artifact
COPY . /artifact
RUN JOBS="$(nproc)" ae/build_all.sh
CMD ["/bin/bash"]
