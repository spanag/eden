# building stuff in docker

DOCKER_TESTING_DIR := $(abspath $(dir $(firstword $(MAKEFILE_LIST) ) ) )
TESTING_DIR ?= $(DOCKER_TESTING_DIR)/../../
PROJ_BASE ?= $(TESTING_DIR)/../

USER_UID_GID = $(shell id -u):$(shell id -g)

EDEN_DOCKER_BUILD_ENV_IMAGE_NAME?=eden-build
EDEN_DOCKER_BUILD_MANYLINUX_IMAGE_NAME?=eden-manylinux
EDEN_DOCKER_TEST_ENV_IMAGE_NAME?=eden-test
EDEN_DOCKER_IMAGE_NAME?=eden-standalone
EDEN_DOCKER_IMAGE_TAG?=latest

EDEN_DOCKER_BUILD_MANYLINUX_IMAGE_FULLNAME?=${EDEN_DOCKER_BUILD_MANYLINUX_IMAGE_NAME}:${EDEN_DOCKER_IMAGE_TAG}
EDEN_DOCKER_BUILD_ENV_IMAGE_FULLNAME?=${EDEN_DOCKER_BUILD_ENV_IMAGE_NAME}:${EDEN_DOCKER_IMAGE_TAG}
EDEN_DOCKER_TEST_ENV_IMAGE_FULLNAME?=${EDEN_DOCKER_TEST_ENV_IMAGE_NAME}:${EDEN_DOCKER_IMAGE_TAG}
EDEN_DOCKER_IMAGE_FULLNAME?=${EDEN_DOCKER_IMAGE_NAME}:${EDEN_DOCKER_IMAGE_TAG}

SUDO_DOCKER:=${DOCKER_TESTING_DIR}/sudo_docker.bash # workaround to omit sudo if not needed to call docker

docker_build_env: ${DOCKER_TESTING_DIR}/build_env.Dockerfile
	${SUDO_DOCKER} build -t ${EDEN_DOCKER_BUILD_ENV_IMAGE_FULLNAME} - < $<
docker_build_manylinux: ${DOCKER_TESTING_DIR}/build_manylinux.Dockerfile
	${SUDO_DOCKER} build --build-arg MANYLINUX_IMAGE=${MANYLINUX_IMAGE} -t ${EDEN_DOCKER_BUILD_MANYLINUX_IMAGE_FULLNAME} - < $<

scipy-notebook-plus: ${DOCKER_TESTING_DIR}/scipy-notebook-plus.Dockerfile
	${SUDO_DOCKER} build -t scidash/scipy-notebook-plus - < $<

neuron-mpi-neuroml: ${DOCKER_TESTING_DIR}/neuron-mpi-neuroml.Dockerfile scipy-notebook-plus
	${SUDO_DOCKER} build -t scidash/neuron-mpi-neuroml - < $< 

docker_test_env: ${DOCKER_TESTING_DIR}/test_env.Dockerfile neuron-mpi-neuroml
	${SUDO_DOCKER} build -t ${EDEN_DOCKER_TEST_ENV_IMAGE_FULLNAME} --file $< ${PROJ_BASE}

build_eden_for_docker: ${DOCKER_TESTING_DIR}/eden_standalone.Dockerfile docker_build_env
	${SUDO_DOCKER} build --file $< -t ${EDEN_DOCKER_IMAGE_FULLNAME} ${PROJ_BASE}

# automatic, containerized testing
# add all containers as dependencies to ensure that at least, they are built, LATER actually run the tests for all images
test: docker_test_env build_eden_for_docker
	rm -rf testing/sandbox/validation_tests
	cp -rT ${PROJ_BASE}testing/validation_tests testing/sandbox/validation_tests
	chmod -R 777 testing/sandbox/validation_tests
	
	${SUDO_DOCKER} run -it --rm --mount type=bind,source=${PROJ_BASE}/testing/sandbox/validation_tests,destination=/home/jovyan/work --group-add 100 --user ${USER_UID_GID} ${EDEN_DOCKER_TEST_ENV_IMAGE_FULLNAME}  python validation_tests.py
	
	rm -rf testing/sandbox/validation_tests
	
	
