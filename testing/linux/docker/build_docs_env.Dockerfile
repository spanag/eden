FROM jupyter/scipy-notebook:2023-08-25
# Build environment for EDEN simulator
MAINTAINER Sotirios Panagiotou <info@sotiriospanagiotou.com>

ENV TEMPREPO /repo_
COPY ./.binder/apt.txt $TEMPREPO/.binder/

USER jovyan
RUN wget -q -N https://www.googleapis.com/download/storage/v1/b/chromium-browser-snapshots/o/Linux_x64%2F1047731%2Fchrome-linux.zip\?generation\=1663284576100523\&alt\=media -O  $HOME/chrome.zip
USER root
RUN chmod 777 -R  $TEMPREPO

# Get the necessary packages to run the docs examples
# perhaps use a ENV PACKAGES variable LATER
# https://www.monolune.com/articles/installing-apt-packages-from-a-requirements-file/
RUN apt-get update \
&& sed 's/#.*//' "$TEMPREPO/.binder/apt.txt" | xargs apt-get install -y \
&& apt-get clean && rm -rf /var/cache/apt/* && rm -rf /var/lib/apt/lists/* && rm -rf /tmp/* 

# TODO on separate apt file i guess...
# from https://github.com/readthedocs/readthedocs-docker-images/blob/ubuntu-22.04-2024.01.29/Dockerfile
# LaTeX -- split to reduce image layer size
RUN apt-get update 
RUN apt-get -y install \
    texlive-fonts-extra
RUN apt-get -y install \
    texlive-lang-english \
    texlive-lang-japanese
RUN apt-get -y install \
    texlive-full

# lmodern: extra fonts
# https://github.com/rtfd/readthedocs.org/issues/5494
#
# xindy: is useful to generate non-ascii indexes
# https://github.com/rtfd/readthedocs.org/issues/4454
#
# fonts-noto-cjk-extra
# fonts-hanazono: chinese fonts
# https://github.com/readthedocs/readthedocs.org/issues/6319
RUN apt-get -y install \
    fonts-symbola \
    lmodern \
    latex-cjk-chinese-arphic-bkai00mp \
    latex-cjk-chinese-arphic-gbsn00lp \
    latex-cjk-chinese-arphic-gkai00mp \
    texlive-fonts-recommended \
    fonts-noto-cjk-extra \
    fonts-hanazono \
    xindy

USER jovyan
# copy chown looks to the host only for username, oh well
COPY --chown=1000 ./.binder/requirements.txt $TEMPREPO/.binder/
COPY --chown=1000 ./docs/requirements.txt $TEMPREPO/docs/
RUN pip install -r  $TEMPREPO/.binder/requirements.txt

# Also apply binder postbuild
COPY --chown=1000 ./.binder/postBuild $TEMPREPO/.binder/
# RUN chmod 777 -R $TEMPREPO

RUN DONT_BUILD_EDEN=1 bash $TEMPREPO/.binder/postBuild
# && rm -rf /var/cache/apt/* && rm -rf /var/lib/apt/lists/* && rm -rf /tmp/*

USER root
# COPY/ADD are being cute and wil spill contents of folders together, always
# the only cure is to copy to each tgt folder separately :c
COPY       eden/ /repo/eden
COPY    testing/ /repo/testing/
COPY thirdparty/ /repo/thirdparty/
COPY    .binder/ /repo/.binder/
COPY    VERSION/ /repo/VERSION
COPY   Makefile/ /repo/Makefile
RUN chmod 777 -R /repo
# or use buildkit or use .dockerignore

USER jovyan
RUN ls -la /repo
RUN ONLY_BUILD_EDEN=1 bash /repo/.binder/postBuild

# RUN python3 -m pip install -r  $TEMPREPO/.binder/requirements.txt
# no files to copy from build context

# TODO
USER root 
RUN apt-get -y install librsvg2-bin
RUN apt-get -y install texlive-latex-recommended texlive-science
USER jovyan
# how about texlive-full ...

CMD ["bash"]
