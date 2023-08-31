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
RUN python3 -m pip install netpyne==1.0.4.2

#TODO move h5 support to neuron-mpi-neuroml, or even notebook-plus
# HDF5 support on JupyterLab
# (Deprecated) Installing extensions with the jupyter labextension install command is now deprecated and will be removed in a future major version of JupyterLab.
# Users should manage prebuilt extensions with package managers like pip and conda, and extension authors are encouraged to distribute their extensions as prebuilt packages
RUN python3 -m pip install hdf5plugin jupyterlab_hdf hdf5widget jupyterlab_h5web[full]==10.0.0 jupyter-server==2.7.2
# RUN jupyter labextension install @jupyterlab/hdf5 && \
    # fix-permissions /etc/jupyter/

USER root
# Get the necessary build & run tools
# More options:
# m4 autoconf for automakeable projects
# ca-certificates for non-canon package repos
# curl cvs svn git for self-hosted repos
# perhaps use a ENV PACKAGES variable LATER
RUN apt-get update && apt-get install -y \
build-essential gcc \
&& apt-get clean && rm -rf /var/cache/apt/* && rm -rf /var/lib/apt/lists/* && rm -rf /tmp/* 


USER $NB_USER
RUN python3 -m pip install -U pip setuptools auditwheel twine

# Now break out in a separate stage to avoid having to rebuild the image whenever sth changes in the source tree, during dev
# Only the final artifacts actually need to be copied to the image, they are on ${EDEN_INSTALL_DIR} 
FROM before-build AS build-from-source
ARG EDEN_CODE_REPO=/repo
ARG EDEN_INSTALL_DIR
USER root
RUN apt-get update && apt-get install -y \
build-essential gcc \
flex=2.6* bison=2:3* xxd \
&& apt-get clean && rm -rf /var/cache/apt/* && rm -rf /var/lib/apt/lists/* && rm -rf /tmp/* 


#TODO copy just what's required for the build (exclude testing scripts that are not required perhaps?)
USER root

# Assume build context is project root
COPY . ${EDEN_CODE_REPO}
WORKDIR ${EDEN_CODE_REPO}

# Add some basic Python integration
ENV OUT_DIR ${EDEN_INSTALL_DIR}
RUN TARGETS="eden hollow_wheel" bash ./testing/linux/docker/build_on_docker.bash

# one more time, for MPI
ENV OUT_DIR ${EDEN_INSTALL_DIR}_MPI
RUN USE_MPI=1 WHEEL_VERSION=$(cat VERSION) bash ./testing/linux/docker/build_on_docker.bash
RUN rm -r ${EDEN_CODE_REPO}

FROM before-build
USER root
ARG EDEN_INSTALL_DIR
COPY --from=build-from-source ${EDEN_INSTALL_DIR} ${EDEN_INSTALL_DIR}
COPY --from=build-from-source ${EDEN_INSTALL_DIR}_MPI ${EDEN_INSTALL_DIR}_MPI
RUN ln -s ${EDEN_INSTALL_DIR}/bin/eden.release.gcc.cpu.x /usr/local/bin/eden
RUN ln -s ${EDEN_INSTALL_DIR}_MPI/bin/eden.release.gcc.cpu.x /usr/local/bin/eden-mpi

USER $NB_USER
# Install the Python package but use the installed binary
RUN python3 -m pip install ${EDEN_INSTALL_DIR}/bin/*.whl

WORKDIR $HOME

# recently in scipy-notebook, permissions of /home/jovyan changed from 775 to 770.
# To run container as not root (ie not using NB_UID which is the recommended way), either --add-group users (another suggested way) or override permissions in the image including x for stat
# RUN chmod a+rx $HOME $HOME/*
# ------------> Ready to run

USER $NB_USER
WORKDIR $WORK_HOME
# CMD ["start-notebook.sh"]
