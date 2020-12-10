FROM gcr.io/finngen-refinery-dev/bioinformatics:0.5

ADD META_ANALYSIS /META_ANALYSIS

RUN apt-get update && apt-get install -y git

# need to upgrade pip to install jass deps
RUN pip3 install --upgrade pip 

RUN pip3 install git+https://gitlab.pasteur.fr/statistical-genetics/JASS_Pre-processing

RUN pip3 install git+https://gitlab.pasteur.fr/statistical-genetics/raiss.git
