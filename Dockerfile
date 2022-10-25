FROM nvidia/cuda:10.2-cudnn8-runtime-centos7 AS guppy
LABEL Name=cawlr Version=0.0.1

RUN yum makecache \
	&& yum -y update ca-certificates \
	&& yum install -y epel-release \
	&& yum clean all \
	&& yum install -y gcc gcc-c++ gcc-gfortran make perl git wget tar \
	zlib-devel ncurses-devel ncurses bzip2 bzip2-devel xz-devel libcurl-devel \
	python3 python3-devel libjpeg-turbo libtiff-devel libjpeg-devel openjpeg2-devel zlib-devel \
	freetype-devel lcms2-devel libwebp-devel tcl-devel tk-devel \
	harfbuzz-devel fribidi-devel libraqm-devel libimagequant-devel libxcb-devel \
	&& wget https://cdn.oxfordnanoportal.com/software/analysis/ont-guppy-6.1.7-1.el7.x86_64.rpm \
	&& yum install -y ont-guppy-6.1.7-1.el7.x86_64.rpm \
	&& rm ont-guppy-6.1.7-1.el7.x86_64.rpm \
	&& yum clean all


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

RUN pip3 install --user matplotlib numpy scikit-learn pyBigWig

WORKDIR /cawlr
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain nightly --profile minimal
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
RUN cargo build --release --bins -Z unstable-options --out-dir /tools

WORKDIR /tools
RUN cp /cawlr/notebooks/*.py . && chmod +x ./*

FROM guppy as dev
COPY --from=builder /tools /tools
ENV PATH="/tools:${PATH}"
ENV PATH="/scripts:${PATH}"

WORKDIR /
CMD ["/bin/bash"]