# neuron-mpi-neuroml
# author Russell Jarvis rjjarvis@asu.edu
# author Rick Gerkin rgerkin@asu.edu

# 2023-08 version
FROM jupyter/scipy-notebook:2023-08-25
USER root
# RUN chown -R $NB_USER $HOME

#Get a whole lot of GNU core development tools
#version control java development, maven
#Libraries required for building MPI from source
#Libraries required for building NEURON from source
RUN apt-get update
RUN apt-get -y install apt-transport-https ca-certificates
RUN apt-get -y install apt-transport-https curl
RUN apt-get install -y wget bzip2 ca-certificates \
    libglib2.0-0 libxext6 libsm6 libxrender1 \
    git gcc g++ build-essential \
    libncurses-dev openmpi-bin openmpi-doc libopenmpi-dev \
    default-jre default-jdk maven emacs \
    libxml2-dev libxslt-dev \
    python3-venv \
    sudo # TODO rm lists?

RUN python3 -m pip install --upgrade pip

# conda is modern enough
# Make sure every Python file belongs to jovyan
RUN find /opt/conda ! -user $NB_USER -print0 | xargs -0 -I {} chown -h $NB_USER {}
# Remove dangling symlinks
RUN find -L /opt/conda -type l -delete
# Make sure every Python file is writable
RUN find /opt/conda ! -writable -print0 | xargs -0 -I {} chmod 744 {}

USER $NB_USER
ENV WORK_HOME /$HOME/work
WORKDIR $WORK_HOME
