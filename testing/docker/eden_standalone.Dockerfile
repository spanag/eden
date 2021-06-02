FROM eden-build:latest

# Standalone build of for EDEN simulator
MAINTAINER Sotirios Panagiotou <info@sotiriospanagiotou.com>

# Assume build context is project root

#TODO copy just what's required for the build (exclude testing scripts that are not required perhaps?)
COPY . /repo

WORKDIR /repo

ENV OUT_DIR /app
RUN bash ./testing/docker/build_on_docker.bash

FROM debian:buster-20191014-slim

# Get the necessary runtime tools
# including GCC with OpenMP
# perhaps use a ENV PACKAGES variable LATER
RUN apt-get update \
    && apt-get install -y \
	build-essential gcc-7 \
    && apt-get clean && rm -rf /var/cache/apt/* && rm -rf /var/lib/apt/lists/* && rm -rf /tmp/* 

WORKDIR /app

COPY --from=0 /app/bin /app/bin

CMD ["bash"]
