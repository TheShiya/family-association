FROM gaow/base-notebook:1.0.0

MAINTAINER Diana Cornejo  <dmc2245@cumc.columbia.edu>

#Compile seqsimlav.2.9.1 from source code and download the C++ library eigen

USER root

RUN cd /tmp && \
    curl -fsSL https://master.dl.sourceforge.net/project/seqsimla/SeqSIMLA2_9_1.tar.gz -o seqsimla2.tar.gz && \
    tar -xvf seqsimla2.tar.gz && \
    mv SeqSIMLA2_9_1/bin/SeqSIMLA* /usr/local/bin && \
    chmod a+x /usr/local/bin/SeqSIMLA && \
    chmod a+x /usr/local/bin convert && \
    cd - && \
    rm -rf /tmp/*

USER jovyan
