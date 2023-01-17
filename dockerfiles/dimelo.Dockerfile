FROM nvidia/cuda:11.8.0-cudnn8-devel-ubuntu20.04 AS guppy
LABEL Name=cawlr Version=0.0.1

RUN apt update
RUN DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y install tzdata

RUN apt install -y wget lsb-release python3-dev python3-pip libbz2-dev liblzma-dev git \
    && export PLATFORM=$(lsb_release -cs) \
    && wget -O- https://cdn.oxfordnanoportal.com/apt/ont-repo.pub | apt-key add - \
    && echo "deb http://cdn.oxfordnanoportal.com/apt ${PLATFORM}-stable non-free" | tee /etc/apt/sources.list.d/nanoporetech.sources.list \
    && apt update \
    && wget -nv https://cdn.oxfordnanoportal.com/software/analysis/ont_guppy_6.1.7-1~focal_amd64.deb \
    && apt-get install -y ./ont_guppy_6.1.7-1~focal_amd64.deb \
    && rm ./ont_guppy_6.1.7-1~focal_amd64.deb

RUN apt-get install -y samtools bedtools

RUN pip3 install ont-pyguppy-client-lib ont-fast5-api megalodon

WORKDIR /rerio
RUN git clone https://github.com/nanoporetech/rerio \
    && rerio/download_model.py rerio/basecall_models/res_dna_r941_min_modbases-all-context_v001 \
    && cp rerio/basecall_models/res_dna_r941_min_modbases-all-context_v001 /opt/ont/guppy/data/ \
    && cp rerio/basecall_models/res_dna_r941_min_modbases-all-context_v001.cfg /opt/ont/guppy/data/ \
    && cp rerio/basecall_models/res_dna_r941_min_modbases-all-context_v001.jsn /opt/ont/guppy/data/ \
    && cp -r rerio/basecall_models/barcoding/* /opt/ont/guppy/data/barcoding/

WORKDIR /winnowmap
RUN wget https://github.com/marbl/Winnowmap/archive/refs/tags/v2.03.tar.gz \
    && tar zxvf v2.03.tar.gz \
    && cd Winnowmap-2.03 \
    && make -j8

WORKDIR /minimap2
RUN git clone -b v2.24 https://github.com/lh3/minimap2.git
RUN cd minimap2 && make

ENV PATH="$PATH:/winnowmap/Winnowmap-2.03/bin"
ENV PATH="$PATH:/minimap2/minimap2"

WORKDIR /
CMD ["/bin/bash"]