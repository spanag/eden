FROM debian:buster-20191014-slim
# Build environment for EDEN simulator
MAINTAINER Sotirios Panagiotou <info@sotiriospanagiotou.com>

# Get the necessary build tools
# perhaps use a ENV PACKAGES variable LATER
RUN apt-get update \
&& apt-get install -y \
build-essential gcc-7=7.4.0* \
flex=2.6.4* bison=2:3.3* \
xxd \
python3=3.7.3* \
&& apt-get clean && rm -rf /var/cache/apt/* && rm -rf /var/lib/apt/lists/* && rm -rf /tmp/* 
# More options:
# m4 autoconf for automakeable projects
# ca-certificates for non-canon package repos
# curl cvs svn git for self-hosted repos
# TODO remove python when an independent testing env is specified

WORKDIR /app
# no files to copy from build context

CMD ["bash"]
