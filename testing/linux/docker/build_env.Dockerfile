FROM debian:trixie-20251208-slim
# Build environment for EDEN simulator
MAINTAINER Sotirios Panagiotou <info@sotiriospanagiotou.com>

# Get the necessary build tools
# perhaps use a ENV PACKAGES variable LATER
# surely there's not much variability in package versions for a distro release, ignore them unless needed for each distro version
RUN apt-get update \
&& apt-get install -y \
build-essential gcc \
flex bison \
xxd \
python3 python3-pip 
# \
#&& apt-get clean && rm -rf /var/cache/apt/* && rm -rf /var/lib/apt/lists/* && rm -rf /tmp/* 
# More options:
# m4 autoconf for automakeable projects
# ca-certificates for non-canon package repos
# curl cvs svn git for self-hosted repos

# Python is included in order to build wheels
RUN apt-get install -y \
	python3-virtualenv python3-setuptools python3-wheel 

# NB: --break-system-packages falls back to the old fragile behaviour, if you don't like it, set up a venv for root i guess https://veronneau.org/python-311-pip-and-breaking-system-packages.html
RUN pip install --break-system-packages auditwheel patchelf

RUN python3 -m auditwheel

# RUN apt-get install patchelf>=0.14
RUN patchelf --version; echo; echo

WORKDIR /app
# no files to copy from build context

CMD ["bash"]
