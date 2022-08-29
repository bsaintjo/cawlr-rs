FROM nvidia/cuda:10.2-cudnn8-devel-centos7 AS guppy
LABEL Name=cawlr Version=0.0.1
RUN yum makecache
RUN yum -y update ca-certificates
RUN yum install -y epel-release
RUN yum install -y gcc gcc-c++ gcc-gfortran make perl git wget tar \
	zlib-devel ncurses-devel ncurses bzip2 bzip2-devel xz-devel libcurl-devel

RUN wget https://cdn.oxfordnanoportal.com/software/analysis/ont-guppy-6.1.7-1.el7.x86_64.rpm
RUN yum install -y ont-guppy-6.1.7-1.el7.x86_64.rpm
RUN rm ont-guppy-6.1.7-1.el7.x86_64.rpm

FROM guppy as base
RUN git clone --recursive -b v0.13.3 https://github.com/jts/nanopolish.git
RUN git clone -b v2.24 https://github.com/lh3/minimap2.git
RUN git clone -b 1.16 https://github.com/samtools/samtools.git

WORKDIR /nanopolish
RUN make all
ENV PATH="/nanopolish:${PATH}"

WORKDIR /minimap2
RUN make
ENV PATH="/minimap2:${PATH}"

WORKDIR /
RUN rm -r samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.16/samtools-1.16.tar.bz2
RUN tar xf samtools-1.16.tar.bz2

WORKDIR /samtools-1.16
RUN ./configure
RUN make
RUN make install

FROM base
WORKDIR /cawlr
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain 1.63.0 --profile minimal
ENV PATH="/root/.cargo/bin:${PATH}"
COPY . .
RUN cargo install --path . --locked

WORKDIR /
CMD ["/bin/bash"]