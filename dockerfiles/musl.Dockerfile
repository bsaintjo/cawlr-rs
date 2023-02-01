FROM rust:1.66.1-alpine

RUN apk update && \
    apk add --no-cache musl-dev make perl gfortran libquadmath openblas-static

# RUN apk update 
# RUN apk add --no-cache musl-dev
# RUN apk add --no-cache make
# RUN apk add --no-cache perl
# RUN apk add --no-cache gfortran
# RUN apk add --no-cache libquadmath
# RUN apk add --no-cache openblas-static

WORKDIR /src