FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    g++ \
    cmake \
    make \
    git \
    pkg-config \
    gdb \
    valgrind \
    clinfo \
    ocl-icd-opencl-dev \
    opencl-c-headers \
    opencl-headers \
    libgsl-dev \
    gnuplot \
    doxygen \
    graphviz \
    check \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace

CMD ["/bin/bash"]
