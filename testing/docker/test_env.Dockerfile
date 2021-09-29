FROM scidash/neuron-mpi-neuroml
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


ARG EDEN_CODE_REPO=/repo
ARG EDEN_INSTALL_DIR=/opt/eden
# Assume build context is project root

#TODO copy just what's required for the build (exclude testing scripts that are not required perhaps?)
COPY . ${EDEN_CODE_REPO}

WORKDIR ${EDEN_CODE_REPO}

# Add some basic Python integration
RUN pip install testing/python_package

ENV OUT_DIR ${EDEN_INSTALL_DIR}
RUN bash ./testing/docker/build_on_docker.bash

# one more time, for MPI
ENV OUT_DIR ${EDEN_INSTALL_DIR}_MPI
RUN USE_MPI=1 bash ./testing/docker/build_on_docker.bash

RUN rm -r ${EDEN_CODE_REPO}

RUN ln -s ${EDEN_INSTALL_DIR}/bin/eden.release.gcc.cpu.x /usr/local/bin/eden
RUN ln -s ${EDEN_INSTALL_DIR}_MPI/bin/eden.release.gcc.cpu.x /usr/local/bin/eden-mpi

WORKDIR $HOME


USER $NB_USER
WORKDIR $WORK_HOME

# CMD ["bash"]
# CMD ["start-notebook.sh"]
