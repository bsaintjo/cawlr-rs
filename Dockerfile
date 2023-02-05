FROM nvidia/cuda:11.8.0-cudnn8-devel-rockylinux8 AS guppy
LABEL Name=cawlr Version=0.0.1

RUN dnf makecache --refresh \
	&& dnf -y update ca-certificates \
	&& dnf install -y dnf-plugins-core \
	&& dnf install -y epel-release \
	&& dnf config-manager -y --set-enabled powertools \
	&& dnf update -y \
	&& dnf install -y gcc gcc-c++ gcc-gfortran make perl git wget tar \
	zlib-devel ncurses-devel ncurses bzip2 bzip2-devel xz-devel libcurl-devel \
	python39 python39-devel \
	libzstd-devel libjpeg-turbo libtiff-devel libjpeg-devel openjpeg2 \
	freetype-devel lcms2 libwebp-devel tcl-devel tk-devel \
	harfbuzz-devel fribidi-devel libraqm-devel libimagequant-devel libxcb-devel \
	libaec mold hdf5 hdf5-devel clang autoconf automake \
	&& wget -nv https://cdn.oxfordnanoportal.com/software/analysis/ont-guppy-6.1.7-1.el8.x86_64.rpm \
	&& dnf install -y ont-guppy-6.1.7-1.el8.x86_64.rpm \
	&& rm ont-guppy-6.1.7-1.el8.x86_64.rpm \
	&& dnf autoremove \
	&& dnf clean all

RUN python3.9 -m pip install --user matplotlib numpy scikit-learn pyBigWig pandas seaborn

FROM guppy as base
RUN git clone --recursive -b v0.13.3 https://github.com/jts/nanopolish.git \
	&& git clone -b v1.1 https://github.com/hasindu2008/f5c.git \
	&& git clone -b v2.24 https://github.com/lh3/minimap2.git \
	&& wget -nv https://github.com/samtools/samtools/releases/download/1.16/samtools-1.16.tar.bz2 \
	&& tar xf samtools-1.16.tar.bz2 \
	&& mv samtools-1.16/ samtools/ \
	&& mkdir /tools/

WORKDIR /nanopolish
RUN make all && cp nanopolish /tools/

WORKDIR /f5c
RUN autoreconf --install && ./scripts/install-hts.sh && ./configure && make cuda=1 zstd=1 && cp ./f5c /tools

WORKDIR /minimap2
RUN make && cp ./minimap2 /tools/

WORKDIR /samtools
RUN ./configure && make && make install && cp ./samtools /tools/

WORKDIR /vbz
RUN wget https://github.com/nanoporetech/vbz_compression/releases/download/v1.0.1/ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz && \
	tar -xzf ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz
ENV HDF5_PLUGIN_PATH=/vbz/ont-vbz-hdf-plugin-1.0.1-Linux/usr/local/hdf5/lib/plugin/

WORKDIR /cawlr
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain 1.66 --profile minimal
ENV PATH="/root/.cargo/bin:${PATH}"

# For caching dependencies
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
RUN mold -run cargo build --release --workspace
RUN cp notebooks/*py /tools \
	&& chmod +x /tools/*

FROM guppy as dev
COPY --from=builder /tools /cawlr/target/release/cawlr \
	/cawlr/target/release/analyze-region-mesmlr-detection-pipeline \
	/cawlr/target/release/analyze-region-pipeline \
	/cawlr/target/release/filter_scores \
	/cawlr/target/release/agg-blocks \
	/cawlr/target/release/train-ctrls-pipeline \
	/cawlr/target/release/filter_detection \
	/tools/
ENV PATH="/tools:${PATH}"
ENV PATH="/root/.local/bin:${PATH}"

WORKDIR /
CMD ["/bin/bash"]