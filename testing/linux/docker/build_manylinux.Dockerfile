ARG MANYLINUX_IMAGE=manylinux:latest

FROM $MANYLINUX_IMAGE
ARG MANYLINUX_IMAGE

# Build environment for EDEN simulator based on manulinux images
MAINTAINER Sotirios Panagiotou <info@sotiriospanagiotou.com>

# Get the necessary build tools
# perhaps use a ENV PACKAGES variable LATER
RUN yum update -y \
&& yum install -y \
vim-common wget m4
#&& apt-get clean && rm -rf /var/cache/apt/* && rm -rf /var/lib/apt/lists/* && rm -rf /tmp/* 
# More options:
# m4 autoconf for automakeable projects
# ca-certificates for non-canon package repos
# curl cvs svn git for self-hosted repos

# Install fresh software from upstream
# flex=2.6.4* bison=2:3.3* \
RUN mkdir -p /fromsrc
WORKDIR /fromsrc
RUN set -e; wget -q https://github.com/westes/flex/releases/download/v2.6.4/flex-2.6.4.tar.gz; \
tar -xf flex-*.tar.gz; 
RUN set -e; cd flex-*/; ./configure; make; make install; flex --version
RUN set -e; wget -q http://ftp.gnu.org/gnu/bison/bison-3.8.2.tar.gz; \
tar -xf bison-*.tar.gz; 
RUN set -e; cd bison-*/; ./configure; make; make install; bison --version

WORKDIR /

# Decide on a python3 for the following
ENV python3=python3.9
RUN mkdir /realpython && ln -sfT "$(which $python3)" /realpython/python3
ENV PATH="/realpython:$PATH"

# Python is included in order to build wheels
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
RUN python3 -m pip install -U pip 
RUN python3 -m pip install virtualenv setuptools wheel auditwheel

RUN python3 -m auditwheel
RUN python3 -m pip install patchelf

RUN patchelf --version

WORKDIR /
# no files to copy from build context

CMD ["bash"]
