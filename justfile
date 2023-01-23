docker:
    docker build --pull -f "Dockerfile" -t bsaintjo/cawlr:full "."
    @echo "Image successfully built"

    docker push bsaintjo/cawlr:full
    @echo "Image successfully pushed"