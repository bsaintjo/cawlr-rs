# Docker

## Available Docker images

* `bsaintjo/cawlr:base`

  * `guppy` v6.1.7
  * `nanopolish` v0.13.3
  * `minimap2` v2.24
  * `samtools` v1.16

* `bsaintjo/cawlr:full`
  * everything in `bsaintjo/cawlr:base`
  * `cawlr` v0.4.0

## Example `docker run` command

```bash
docker run --rm \
    -it --gpus all \
    -v "$(pwd)"/path/to/fast5:/fast5 \
    bsaintjo/cawlr:full /bin/bash
```

* The `--gpus all` parameter allows for GPU-accelerated basecalling within the Docker container.
