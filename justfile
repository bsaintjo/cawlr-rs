docker:
    docker build --pull -f "Dockerfile" -t brookslab/cawlr:latest "."
    @echo "Image successfully built"

    docker push brookslab/cawlr:latest
    @echo "Image successfully pushed"

musl:
    docker build --pull -f "dockerfiles/musl.Dockerfile" -t rust-musl-builder:latest "."
    docker run \
        --rm -it \
        -u "$(id -u)":"$(id -g)" \
        -e RUSTFLAGS="-C target-feature=+crt-static" \
        -v "$(pwd)":/src rust-musl-builder:latest \
        cargo build --release --workspace --target x86_64-unknown-linux-musl