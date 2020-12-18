FROM ubuntu:20.04

ADD META_ANALYSIS /META_ANALYSIS

RUN apt-get update && \
    apt-get install -y git && \
    apt-get install python3-pip --yes

# need to upgrade pip to install jass deps
RUN pip3 install --upgrade pip 

RUN pip3 install git+https://gitlab.pasteur.fr/statistical-genetics/JASS_Pre-processing

RUN pip3 install git+https://gitlab.pasteur.fr/statistical-genetics/raiss.git
