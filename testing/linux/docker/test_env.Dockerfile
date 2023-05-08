# So all ARG's that have beed defined in a stage of a multistage Dockerfile have to be
# mentioned by name over and over to exist in the following stages.
# But popular ARG's like NB_USER are always visible even for further stages,
# perhaps because they were stated in the a base image before this Dockerfile.
# Life is often unfair like that.
# See https://github.com/moby/moby/issues/34129#issuecomment-564541933 https://ryandaniels.ca/blog/docker-dockerfile-arg-from-arg-trouble/
ARG EDEN_INSTALL_DIR=/opt/eden

FROM scidash/neuron-mpi-neuroml AS before-build 
# Test & verification environment for EDEN simulator, including GCC, Python3, NEURON, jNeuroML, the works.
# TODO make an independent Dockerfile from ground up, to ensure compatible tool versions.
MAINTAINER Sotirios Panagiotou <info@sotiriospanagiotou.com>

USER $NB_USER

# Get some auxiliar packages, for network generation on the testing environment
# NetPyNE
RUN python3 -m pip install netpyne==1.0.0.2

# HDF5 support on JupyterLab
RUN python3 -m pip install hdf5plugin jupyterlab_hdf hdf5widget
RUN jupyter labextension install @jupyterlab/hdf5 && \
    fix-permissions /etc/jupyter/

USER root

# Get the necessary build tools
# perhaps use a ENV PACKAGES variable LATER
RUN apt-get update \
&& apt-get install -y \
build-essential gcc-7 \
flex=2.6* bison=2:3* \
xxd \
&& apt-get clean && rm -rf /var/cache/apt/* && rm -rf /var/lib/apt/lists/* && rm -rf /tmp/* 
# More options:
# m4 autoconf for automakeable projects
# ca-certificates for non-canon package repos
# curl cvs svn git for self-hosted reposd

# Get the necessary runtime tools
# including GCC with OpenMP
# perhaps use a ENV PACKAGES variable LATER
RUN apt-get update \
    && apt-get install -y \
	build-essential gcc-7 \
    && apt-get clean && rm -rf /var/cache/apt/* && rm -rf /var/lib/apt/lists/* && rm -rf /tmp/* 

USER $NB_USER
RUN python3 -m pip install -U pip setuptools auditwheel twine

# Assume build context is project root

# Now break out in a separate stage to avoid having to rebuild the image whenever sth changes in the source tree, during dev
# Only the final artifacts actually need to be copied to the image, they are on ${EDEN_INSTALL_DIR} 
FROM before-build AS build-from-source
ARG EDEN_CODE_REPO=/repo
ARG EDEN_INSTALL_DIR

#TODO copy just what's required for the build (exclude testing scripts that are not required perhaps?)
USER root
# RUN echo $NB_USER
COPY . ${EDEN_CODE_REPO}

WORKDIR ${EDEN_CODE_REPO}

# Add some basic Python integration
# RUN pip install testing/python_package

ENV OUT_DIR ${EDEN_INSTALL_DIR}
RUN TARGETS="eden hollow_wheel" bash ./testing/linux/docker/build_on_docker.bash

# one more time, for MPI
ENV OUT_DIR ${EDEN_INSTALL_DIR}_MPI
RUN USE_MPI=1 WHEEL_VERSION=$(cat VERSION) bash ./testing/linux/docker/build_on_docker.bash

FROM before-build
ARG EDEN_INSTALL_DIR
COPY --from=build-from-source ${EDEN_INSTALL_DIR} ${EDEN_INSTALL_DIR}
COPY --from=build-from-source ${EDEN_INSTALL_DIR}_MPI ${EDEN_INSTALL_DIR}_MPI
# RUN rm -r ${EDEN_CODE_REPO}

USER $NB_USER
# Install the Python package but use the installed binary
RUN python3 -m pip install ${EDEN_INSTALL_DIR}/bin/*.whl

USER root

RUN ln -s ${EDEN_INSTALL_DIR}/bin/eden.release.gcc.cpu.x /usr/local/bin/eden
RUN ln -s ${EDEN_INSTALL_DIR}_MPI/bin/eden.release.gcc.cpu.x /usr/local/bin/eden-mpi

WORKDIR $HOME


# ------------> Ready to run

USER $NB_USER
WORKDIR $WORK_HOME

# CMD ["bash"]
# CMD ["start-notebook.sh"]
