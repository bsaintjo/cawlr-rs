
FROM nvidia/cuda:11.8.0-cudnn8-devel-ubuntu20.04 AS guppy
LABEL Name=cawlr Version=0.0.1

RUN apt update
RUN DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y install tzdata
RUN apt install -y wget lsb-release python3-dev python3-pip gfortran sqlite3 apt-utils \
    git autoconf libhdf5-dev libzstd-dev curl \
    && export PLATFORM=$(lsb_release -cs) \
    && wget -O- https://cdn.oxfordnanoportal.com/apt/ont-repo.pub | apt-key add - \
    && echo "deb http://cdn.oxfordnanoportal.com/apt ${PLATFORM}-stable non-free" | tee /etc/apt/sources.list.d/nanoporetech.sources.list \
    && apt update \
    && wget -nv https://cdn.oxfordnanoportal.com/software/analysis/ont_guppy_6.1.7-1~focal_amd64.deb \
    && apt-get install -y ./ont_guppy_6.1.7-1~focal_amd64.deb \
    && rm ./ont_guppy_6.1.7-1~focal_amd64.deb
RUN python3 -m pip install --user matplotlib numpy scikit-learn pyBigWig pandas seaborn

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
RUN ./configure --without-curses --disable-bz2 --disable-lzma && make && make install && cp ./samtools /tools/

WORKDIR /cawlr
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain 1.66 --profile minimal
ENV PATH="/root/.cargo/bin:${PATH}"

# For caching dependencies
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain 1.66 --profile minimal
ENV PATH="/root/.cargo/bin:${PATH}"
RUN cargo install cargo-chef

FROM base AS planner
WORKDIR /cawlr
COPY . .
RUN cargo chef prepare --recipe-path recipe.json

FROM base AS builder
WORKDIR /cawlr
COPY --from=planner /cawlr/recipe.json recipe.json
RUN wget https://github.com/rui314/mold/releases/download/v1.7.1/mold-1.7.1-x86_64-linux.tar.gz \
    && tar xf mold-1.7.1-x86_64-linux.tar.gz
RUN cargo chef cook --release --recipe-path recipe.json
COPY . .
RUN mold-1.7.1-x86_64-linux/bin/mold -run cargo build --release --bins --all-features
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