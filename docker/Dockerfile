FROM ubuntu:20.04

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

ADD META_ANALYSIS /META_ANALYSIS

ENV TZ=Europe/Helsinki
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && \
    apt-get install curl make gcc g++ cmake gfortran libreadline-dev \
    libz-dev libzstd-dev libbz2-dev liblzma-dev libpcre3-dev libssl-dev libcurl4-openssl-dev \
    libopenblas-dev default-jre unzip libboost-all-dev \
    libpng-dev libcairo2-dev tabix --yes && \
    apt-get clean &&

RUN curl -O https://cloud.r-project.org/src/base/R-3/R-3.6.1.tar.gz && \
    tar xvzf R-3.6.1.tar.gz && \
    cd R-3.6.1 && \
    ./configure --with-x=no --enable-R-shlib=yes --with-blas="-lopenblas" && \
    make && \
    mkdir -p /usr/local/lib/R/lib && \
    make install && \
    cd .. && \
    rm -rf R-3.6.1*

RUN apt-get update && \
    apt-get install -y git && \
    apt-get install python3-pip --yes

# need to upgrade pip to install jass deps
RUN pip3 install --upgrade pip 

RUN pip3 install git+https://gitlab.pasteur.fr/statistical-genetics/JASS_Pre-processing

RUN pip3 install git+https://gitlab.pasteur.fr/statistical-genetics/raiss.git