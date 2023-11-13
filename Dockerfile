FROM ubuntu:20.04

ENV DEBIAN_FRONTEND noninteractive

# Instal basic utilities
RUN apt-get update
RUN apt-get install -y --no-install-recommends git wget unzip bzip2 build-essential ca-certificates pip gnupg default-jdk
RUN apt-get update

## Add the antiSMASH debian repository
RUN apt-get install -y apt-transport-https
RUN wget http://dl.secondarymetabolites.org/antismash-stretch.list -O /etc/apt/sources.list.d/antismash.list
RUN wget -q -O- http://dl.secondarymetabolites.org/antismash.asc | apt-key add -
RUN apt-get update

## Install the binaries
RUN apt-get install -y hmmer2 hmmer diamond-aligner fasttree prodigal ncbi-blast+ muscle glimmerhmm meme-suite
RUN pip install Jinja2==3.0.1

## Install antiSMASH
RUN wget https://dl.secondarymetabolites.org/releases/6.1.1/antismash-6.1.1.tar.gz
RUN tar -zxf antismash-6.1.1.tar.gz
RUN pip install ./antismash-6.1.1
COPY proteins.dmnd /usr/local/lib/python3.8/dist-packages/antismash/databases/clusterblast/proteins.dmnd
RUN download-antismash-databases
RUN antismash --check-prereqs

# Install FunARTS

## Add related files
ADD . /funarts/
RUN mkdir /custom_files/
RUN mkdir /custom_reference/
WORKDIR /funarts/

## Install pip depencies
RUN pip install -r requirements.txt

## Cloning BiG-SCAPE
#RUN git clone https://github.com/medema-group/BiG-SCAPE.git
##https://github.com/medema-group/BiG-SCAPE/archive/refs/tags/v1.1.5.tar.gz
## Install bigscape
#RUN cd BiG-SCAPE \
# && wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz \
# && gunzip Pfam-A.hmm.gz \
# && hmmpress Pfam-A.hmm \
# && chmod +x /funarts/BiG-SCAPE/*py \
# && chmod a+w /funarts/BiG-SCAPE/domains_color_file.tsv \
# && chmod a+w /funarts/BiG-SCAPE/Annotated_MIBiG_reference/

RUN cd BiG-SCAPE \
 && chmod +x /funarts/BiG-SCAPE/*py \
 && chmod a+w /funarts/BiG-SCAPE/domains_color_file.tsv \
 && chmod a+w /funarts/BiG-SCAPE/Annotated_MIBiG_reference/

## Remove unnecessery files
RUN rm proteins.dmnd \
 && rm /antismash-6.1.1.tar.gz \
 && rm -rf /var/lib/apt/lists/*

ENTRYPOINT ["python3", "funartspipeline.py"]

