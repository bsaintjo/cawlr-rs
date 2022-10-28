FROM nvidia/cuda:11.8.0-cudnn8-runtime-rockylinux8 AS guppy
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
	libjpeg-turbo libtiff-devel libjpeg-devel openjpeg2 \
	zlib-devel freetype-devel lcms2 libwebp-devel tcl-devel tk-devel \
	harfbuzz-devel fribidi-devel libraqm-devel libimagequant-devel libxcb-devel \
	libaec mold hdf5 hdf5-devel clang \
	&& wget -nv https://cdn.oxfordnanoportal.com/software/analysis/ont-guppy-6.1.7-1.el8.x86_64.rpm \
	&& dnf install -y ont-guppy-6.1.7-1.el8.x86_64.rpm \
	&& rm ont-guppy-6.1.7-1.el8.x86_64.rpm \
	&& dnf autoremove \
	&& dnf clean all

RUN python3.9 -m pip install --user matplotlib numpy scikit-learn pyBigWig pandas seaborn

FROM guppy as base
RUN git clone --recursive -b v0.13.3 https://github.com/jts/nanopolish.git \
	&& git clone -b v2.24 https://github.com/lh3/minimap2.git \
	&& wget -nv https://github.com/samtools/samtools/releases/download/1.16/samtools-1.16.tar.bz2 \
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
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain nightly --profile minimal
ENV PATH="/root/.cargo/bin:${PATH}"

# For caching depedencies
RUN mold -run cargo install cargo-chef
FROM base AS planner
WORKDIR /cawlr
COPY . .
RUN mold -run cargo chef prepare --recipe-path recipe.json

FROM base AS builder
WORKDIR /cawlr
COPY --from=planner /cawlr/recipe.json recipe.json
RUN mold -run cargo chef cook --release --recipe-path recipe.json
COPY ./*[^py] ./
RUN mold -run cargo build --release --all-features --bins -Z unstable-options --out-dir /tools
COPY ./notebooks/*py ./notebooks/
RUN cp notebooks/*py /tools \
	&& chmod +x /tools/*

FROM guppy as dev
COPY --from=builder /tools /tools
ENV PATH="/tools:${PATH}"
ENV PATH="/root/.local/bin:${PATH}"

WORKDIR /
CMD ["/bin/bash"]