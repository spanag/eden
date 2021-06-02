
# Assume runnning from Parallel_HH/ directory, on Docker container

# Default parameters

# Make Targets to build
: "${TARGETS:=eden}"

# Directories for Artifacts
: "${OUT_DIR:=$(mktemp -d)}"
mkdir -p ${OUT_DIR}/bin
mkdir -p ${OUT_DIR}/obj

TOOLCHAIN=gcc OUT_DIR=${OUT_DIR} BUILD=release make -j$(nproc) ${TARGETS}

if [ $? -ne 0 ]; then
	echo "Dockerized build with GCC failed !"
	exit 1
fi
