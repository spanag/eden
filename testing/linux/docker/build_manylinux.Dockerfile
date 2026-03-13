ARG MANYLINUX_IMAGE=manylinux:latest

FROM $MANYLINUX_IMAGE
ARG MANYLINUX_IMAGE

# Build environment for EDEN simulator based on manulinux images
MAINTAINER Sotirios Panagiotou <info@sotiriospanagiotou.com>


# for manylinux2014: Centos 7 is EOL and is no longer available from the usual mirrors, so switch to https://vault.centos.org
# see https://github.com/pypa/manylinux/pull/1628/commits/7beb9ae220bcf3da425d323817709c1a1e2bd35d https://github.com/pypa/manylinux/issues/1641
RUN if [[ "${MANYLINUX_IMAGE}" == *"manylinux2014"* ]]; then \
	sed -i 's/enabled=1/enabled=0/g' /etc/yum/pluginconf.d/fastestmirror.conf; \
	sed -i 's/^mirrorlist/#mirrorlist/g' /etc/yum.repos.d/*.repo; \
	sed -i 's;^.*baseurl=http://mirror;baseurl=https://vault;g' /etc/yum.repos.d/*.repo; \
	if [ "${AUDITWHEEL_ARCH}" == "aarch64" ] || [ "${AUDITWHEEL_ARCH}" == "ppc64le" ]; then \
		sed -i 's;/centos/7/;/altarch/7/;g' /etc/yum.repos.d/*.repo ;\
	fi;\
fi



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
