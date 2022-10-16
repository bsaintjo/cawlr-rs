FROM nvidia/cuda:10.2-cudnn8-devel-centos7 AS guppy
LABEL Name=cawlr Version=0.0.1

RUN yum makecache \
	&& yum -y update ca-certificates \
	&& yum install -y epel-release \
	&& yum install -y gcc gcc-c++ gcc-gfortran make perl git wget tar \
	zlib-devel ncurses-devel ncurses bzip2 bzip2-devel xz-devel libcurl-devel \
	&& wget https://cdn.oxfordnanoportal.com/software/analysis/ont-guppy-6.1.7-1.el7.x86_64.rpm \
	&& yum install -y ont-guppy-6.1.7-1.el7.x86_64.rpm \
	&& rm ont-guppy-6.1.7-1.el7.x86_64.rpm

FROM guppy as base
RUN git clone --recursive -b v0.13.3 https://github.com/jts/nanopolish.git \
	&& git clone -b v2.24 https://github.com/lh3/minimap2.git \
	&& wget https://github.com/samtools/samtools/releases/download/1.16/samtools-1.16.tar.bz2 \
	&& tar xf samtools-1.16.tar.bz2 \
	&& mv samtools-1.16/ samtools/ \
	&& mkdir /tools/

WORKDIR /nanopolish
RUN make all && cp nanopolish /tools/
ENV PATH="/nanopolish:${PATH}"

WORKDIR /minimap2
RUN make && cp ./minimap2 /tools/

WORKDIR /samtools
RUN ./configure && make && make install && cp ./samtools /tools/

WORKDIR /cawlr
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain stable --profile minimal
ENV PATH="/root/.cargo/bin:${PATH}"

# For caching depedencies
RUN cargo install cargo-chef
FROM base AS planner
WORKDIR /cawlr
COPY . .
RUN cargo chef prepare --recipe-path recipe.json

FROM base AS builder
WORKDIR /cawlr
COPY --from=planner /cawlr/recipe.json recipe.json
RUN cargo chef cook --release --recipe-path recipe.json
COPY . .
RUN cargo build --release --bins \
	&& cp ./target/release/cawlr /tools/ \
	&& cp ./target/release/convert-detection /tools/ \
	&& cp ./target/release/analyze-region-pipeline /tools/ \
	&& cp ./target/release/overlap-bed /tools/ \
	&& cp ./target/release/agg-blocks /tools/

FROM guppy as dev
COPY --from=base /tools/ /tools/
ENV PATH="/tools:${PATH}"

WORKDIR /
CMD ["/bin/bash"]