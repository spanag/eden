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
RUN apt-get update && apt-get install -y wget bzip2 ca-certificates default-jre default-jdk maven automake libtool  \
                       wget python3 libpython3-dev libncurses5-dev libreadline-dev libgsl0-dev cython3 \
                       python3-pip python3-numpy python3-scipy python3-matplotlib python3-mock \
                       ipython3 python3-docutils \
                       python3-mpi4py cmake ssh


#Do the rest of the build  as user:
#This will create a more familiar environment to continue developing in.
#with less of a need to chown and chmod everything done as root at dockerbuild completion

USER jovyan
# Use numpy 1.12.1 until quantities is compatible with 1.13.
# RUN conda install -y scipy numpy==1.12.1 matplotlib
RUN sudo chown -R jovyan /home/jovyan
ENV HOME /home/jovyan
ENV PATH /opt/conda/bin:/opt/conda/bin/conda:/opt/conda/bin/python:$PATH

#Test matplotlib
RUN python -c "import matplotlib"
#Install General MPI, such that mpi4py can later bind with it.

WORKDIR $HOME

ARG NEURON_VERSION=7.7
ARG NEURON_VERSIONED=nrn-${NEURON_VERSION}

# could make stderr more quiet or sth, to suppress red text
RUN \
  wget http://www.neuron.yale.edu/ftp/neuron/versions/v${NEURON_VERSION}/${NEURON_VERSIONED}.tar.gz && \
  tar -xzf ${NEURON_VERSIONED}.tar.gz && \
  rm ${NEURON_VERSIONED}.tar.gz

WORKDIR $HOME/${NEURON_VERSIONED}
ENV PATH /usr/bin/python3/python:/opt/conda/bin:/opt/conda/bin/conda:/opt/conda/bin/python:$PATH
RUN ./configure --prefix=`pwd` --with-paranrn --without-iv --with-nrnpython=/opt/conda/bin/python3
RUN sudo make all && \
   make install
   
RUN make all
RUN make install

WORKDIR src/nrnpython
RUN python setup.py install
ENV NEURON_HOME $HOME/${NEURON_VERSIONED}/x86_64
ENV PATH $NEURON_HOME/bin:$PATH


# Install NeuroML tools
USER root
# required for libNeuroML
RUN apt-get install -y python-lxml python3-pip python-dev python3-setuptools

USER jovyan
RUN sudo chown -R jovyan /home/jovyan
WORKDIR $HOME

RUN pip3 install --upgrade pip
RUN pip3 install --upgrade virtualenv

# to get from pypi (not recommended)
# RUN pip install pyNeuroML==0.3.11

# Get the original source code
# Could also download just the branch to be used,
# in case no other branch is needed elsewhere.
# But picking the branch is realted to the specific commit, so it's a bit laborious and error-friendly.
RUN git clone https://github.com/NeuralEnsemble/libNeuroML.git
RUN git clone https://github.com/NeuroML/pyNeuroML.git

# NOTE: building jNeuroML fails, at building the 'injecting plugin'
# so we actually use pyNeuroML's bundled jar instead.
# RUN git clone https://github.com/NeuroML/jNeuroML
# WORKDIR jNeuroML
# RUN git checkout development
# RUN ls -ltr *.py
# RUN python getNeuroML.py

# Build from source code

WORKDIR $HOME/libNeuroML
# 2020-05 version
RUN git checkout 2807173626d3d00e6ffadd728a43465603a54431
RUN pip3 install . -r requirements.txt

WORKDIR $HOME/pyNeuroML
# 2020-05 version
RUN git checkout d34dd127c3875c89700678ddbc2219d35b4f0e7e
RUN python3 -m pip install .


# Check if packages are imported OK
WORKDIR $WORK_HOME

RUN python3 -c "import neuroml"
RUN python3 -c "import neuroml; from pyneuroml import pynml"

RUN ExampleDir=$(mktemp -d); cp -r ~/pyNeuroML/examples/ $ExampleDir; cd $ExampleDir/examples; python run_jneuroml_plot_matplotlib.py -nogui <&-

# TODO more testing on NeuroML examples, to verify it's working properly (like with Poisson sources and such)

#onward
WORKDIR $WORK_HOME
