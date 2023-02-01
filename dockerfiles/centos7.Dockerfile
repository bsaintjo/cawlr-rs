FROM centos:7
LABEL Name=cawlr-centos7 Version=0.0.1
RUN yum install -y epel-release
RUN yum install -y gcc gcc-gfortran make perl git
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain nightly --profile minimal
ENV PATH="/root/.cargo/bin:${PATH}"
WORKDIR /cawlr
COPY . .
RUN cargo build --release --locked
RUN strip "target/release/cawlr"
RUN cargo test --bin cawlr --locked
RUN cargo install --path . --locked
RUN cawlr -h
CMD echo "777"