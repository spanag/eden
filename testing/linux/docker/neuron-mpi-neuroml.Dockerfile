#author russell jarvis rjjarvis@asu.edu
#NEURON Dockerfile
#Set the base image to Ubuntu
FROM scidash/scipy-notebook-plus

#Get a whole lot of GNU core development tools
#version control java development, maven
#Libraries required for building MPI from source
#Libraries required for building NEURON from source
#Also DO this part as root.
USER root
RUN apt-get update && apt-get install -y \
    wget bzip2 ca-certificates default-jre default-jdk maven automake libtool cmake ssh\
    libncurses5-dev libreadline-dev libgsl0-dev \
    cython3 libpython3-dev python3-mock ipython3 python3-docutils \
    python3-mpi4py

#Do the rest of the build as user:
USER jovyan

# Test matplotlib, just check if it is imported ok
RUN python -c "import matplotlib"
#Install General MPI, such that mpi4py can later bind with it.
# ---> MPI is actually already installed in previous image

# ---> Install NEURON
WORKDIR $HOME

ARG NEURON_VERSION=8.2.2
RUN python3 -m pip install neuron==${NEURON_VERSION}
RUN which nrnivmodl nrniv; echo $PATH
# NEURON_HOME is needed by NML tooling
# but ENV can't be determined from inside the container being built https://github.com/moby/moby/issues/29110
# ENV NEURON_HOME $(realpath $(dirname $(which nrniv))/..)
ENV NEURON_HOME=/opt/conda

# ---> Install NeuroML tools
USER root
# required for h5 ?
RUN sudo apt-get install -y libhdf5-dev
# required for libNeuroML ?
# RUN apt-get install -y python3-lxml python3-pip python3-setuptools

USER jovyan
RUN python3 -m pip install tables>=3.3.0

# ---> Install NeuroML tool reqs
USER jovyan

RUN python3 -m pip install libneuroml==0.5.5 pyneuroml==1.1.7
# or ...Get the original source code
# WORKDIR $HOME
# Could also download just the branch to be used,
# in case no other branch is needed elsewhere.
# But picking the branch is related to the specific commit, so it's a bit laborious and error-friendly.
# RUN git clone https://github.com/NeuralEnsemble/libNeuroML.git

# # NOTE: building jNeuroML fails, at building the 'injecting plugin'
# # so we actually use pyNeuroML's bundled jar instead.
# # RUN git clone https://github.com/NeuroML/jNeuroML
# # WORKDIR jNeuroML
# # RUN git checkout development; ls -ltr *.py
# # RUN python getNeuroML.py

# Build from source code
# WORKDIR $HOME/libNeuroML
# 2021-03 version
# RUN git checkout 632c1bce797d44308d5ec8246c0aac360c862f1a
# RUN python3 -m pip install . -r requirements.txt

# WORKDIR $HOME/pyNeuroML
# 2021-03 version
# RUN git checkout 9e070467498c57d4244d44f9996bb8e0eecc5dc3
# RUN python3 -m pip install .RUN 

# Check if packages work OK
RUN python3 -c "import neuroml; from pyneuroml import pynml"
# run an example as well
WORKDIR $HOME
RUN git clone --depth=1 https://github.com/NeuroML/pyNeuroML.git
WORKDIR $WORK_HOME
RUN ExampleDir=$(mktemp -d); cp -r ~/pyNeuroML/examples/ $ExampleDir; cd $ExampleDir/examples; python run_jneuroml_plot_matplotlib.py -nogui <&-

# TODO more testing on NeuroML examples, to verify it's working properly (like with Poisson sources and such)

#onward
WORKDIR $WORK_HOME
