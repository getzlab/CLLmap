FROM ubuntu:16.04
MAINTAINER Binyamin A. Knisbacher

# update package manager & build essentials
RUN apt-get update \
	&& apt-get install --yes build-essential wget libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev \
	automake autoconf

FROM gcc:5.1 

ENV SAMTOOLS_VERSION 1.9

# install samtools
RUN cd /tmp \
    && wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar xvjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && cd samtools-${SAMTOOLS_VERSION} \
    && ./configure --prefix=/usr/local \
    && make \
    && make install \
    && cd / && rm -rf /tmp/samtools-${SAMTOOLS_VERSION} \
    && cd /opt && \
    git clone --recursive https://github.com/jwalabroad/VariantBam.git && \
    cd VariantBam/SeqLib/htslib && \
    autoheader && \
    autoconf && \
    ./configure --enable-libcurl && \
    cd ../../  && \
    ./configure LDFLAGS="-lcurl -lcrypto" && \
    make && \
    make install

ENV PATH /opt/VariantBam/src/variant:$PATH

