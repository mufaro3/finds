FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    gfortran \
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
    xvfb \
    python3 \
    python3.12-venv \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# install FMM3D
COPY lib/FMM3D workspace/lib/FMM3D
RUN cd workspace/lib/FMM3D && make lib

ARG USER_UID=1000
ARG USER_GID=1000

RUN if getent group ${USER_GID}; then \
        group_name=$(getent group ${USER_GID} | cut -d: -f1); \
        groupmod -n dev ${group_name}; \
    else \
        groupadd --gid ${USER_GID} dev; \
    fi \
    && if id -u ${USER_UID} >/dev/null 2>&1; then \
        usermod -l dev -d /home/dev -m $(id -un ${USER_UID}); \
    else \
        useradd --uid ${USER_UID} \
                --gid dev \
                --create-home \
                --shell /bin/bash \
                dev; \
    fi

WORKDIR /workspace

# install requirements
COPY requirements.txt .
RUN pip3 install --break-system-packages --no-cache-dir -r requirements.txt

CMD ["/bin/bash"]
