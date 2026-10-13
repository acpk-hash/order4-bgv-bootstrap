# Asiacrypt 2026 artifact: Order-Four Symmetry in BGV Bootstrapping.
# Toolchain of the measurement host: Ubuntu 22.04, GCC 11.4, NTL 11.5.1
# (default configuration, from 2_baselines/), GMP 6.2.1. Apart from the Ubuntu
# packages, the build uses only files in this repository.
#
#   docker build -t order4-bgv .
#   docker run --rm -it order4-bgv                         # shell in /artifact
#   docker run --rm -v "$PWD/out:/artifact/out" order4-bgv 3_experiments/verify.sh
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential cmake git ca-certificates m4 patch \
        libgmp-dev python3 python3-matplotlib time procps \
    && rm -rf /var/lib/apt/lists/*

# NTL in its default configuration. NTL_ENABLE_AVX_FFT must stay off: it
# shrinks NTL's primes from 60 to 50 bits, which changes the modulus chain
# HElib builds (1066 -> 1086 bits on m = 2^16). check_ntl.cpp enforces this.
COPY 2_baselines/ntl-11.5.1.tar.gz 3_experiments/check_ntl.cpp /tmp/
RUN cd /tmp && tar xzf ntl-11.5.1.tar.gz && cd ntl-11.5.1/src \
    && ./configure SHARED=on && make -j"$(nproc)" && make install \
    && cd /tmp && rm -rf ntl-11.5.1 ntl-11.5.1.tar.gz && ldconfig \
    && g++ /tmp/check_ntl.cpp -o /tmp/check_ntl -lntl -lgmp -lpthread && /tmp/check_ntl

WORKDIR /artifact
COPY . /artifact
RUN JOBS="$(nproc)" 3_experiments/build_all.sh
CMD ["/bin/bash"]
