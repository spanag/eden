#!/bin/bash
set -e
# Assume runnning from Parallel_HH/ directory, on Docker container

# Default parameters
# https://stackoverflow.com/questions/28085062/assigning-default-values-to-shell-variables-with-a-single-command-in-bash
# Make Targets to build
: "${TARGETS:=eden}"

# Directories for Artifacts
: "${OUT_DIR:=$(mktemp -d)}"
mkdir -p ${OUT_DIR}/bin
mkdir -p ${OUT_DIR}/obj

procs=$(nproc)
# cut build parallelism for less capable machines
if [ $(awk '/^MemAvailable:/ { print $2; }' /proc/meminfo) -lt 1000000 ]; then
	procs=1
fi

TOOLCHAIN=gcc OUT_DIR=${OUT_DIR} BUILD=release make -j${procs} ${TARGETS}

if [ $? -ne 0 ]; then
	echo "Dockerized build with GCC failed !"
	exit 1
fi
