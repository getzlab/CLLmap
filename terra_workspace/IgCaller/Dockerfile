FROM ubuntu:16.04
MAINTAINER Binyamin A. Knisbacher

RUN apt-get update && apt-get install -y software-properties-common && \
    apt-get update && apt-get install -y \
        build-essential \
        cmake \
        curl \
        git \
        libboost-all-dev \
        libbz2-dev \
        libcurl3-dev \
        libhdf5-serial-dev \
        liblzma-dev \
        libncurses5-dev \
        libssl-dev \
        python3 \
        python3-pip \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/

# htslib
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
    tar -xf htslib-1.9.tar.bz2 && rm htslib-1.9.tar.bz2 && cd htslib-1.9 && \
    ./configure --enable-libcurl --enable-s3 --enable-plugins --enable-gcs && \
    make && make install && make clean

# samtools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar -xf samtools-1.9.tar.bz2 && rm samtools-1.9.tar.bz2 && cd samtools-1.9 && \
    ./configure --with-htslib=/opt/htslib-1.9 && make && make install && make clean

# python modules
RUN pip3 install --upgrade pip setuptools
RUN pip3 install statistics regex
RUN pip3 install 'argparse==1.1' 'numpy==1.16.3' 'scipy==v1.3.0'

# IgCaller - including creating IgCaller symlink to executable
RUN mkdir /src && cd /src && git clone https://github.com/ferrannadeu/IgCaller.git && \
	/bin/bash -c 'echo -e "#!/usr/bin/env python3" | cat - /src/IgCaller/IgCaller_v1.1.py > /src/IgCaller/IgCaller' && \
	/bin/bash -c 'chmod +x /src/IgCaller/IgCaller'
ENV PATH="/src/IgCaller:${PATH}"
