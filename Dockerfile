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
    libjson-c-dev \
    libhdf5-dev \
    libncurses-dev \
    doxygen \
    graphviz \
    check \
    ffmpeg \
    texlive-latex-base \
    texlive-latex-extra \
    texlive-latex-recommended \
    texlive-fonts-extra \
    texlive-bibtex-extra \
    biber \
    latexmk \
    python3 \
    python3.12-venv \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

USER 1000

WORKDIR /workspace

COPY requirements.txt .

RUN pip3 install --break-system-packages --no-cache-dir -r requirements.txt

CMD ["/bin/bash"]
