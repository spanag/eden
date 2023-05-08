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
python3=3.7.3* python3-pip 
# \
#&& apt-get clean && rm -rf /var/cache/apt/* && rm -rf /var/lib/apt/lists/* && rm -rf /tmp/* 
# More options:
# m4 autoconf for automakeable projects
# ca-certificates for non-canon package repos
# curl cvs svn git for self-hosted repos

# Python is included in order to build wheels
RUN python3 -m pip install -U pip
RUN python3 -m pip install virtualenv setuptools wheel auditwheel

RUN python3 -m auditwheel
RUN python3 -m pip install patchelf

# RUN apt-get install patchelf>=0.14
RUN patchelf --version; echo; echo

WORKDIR /app
# no files to copy from build context

CMD ["bash"]
