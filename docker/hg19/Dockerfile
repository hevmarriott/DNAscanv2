############################################################
# Dockerfile to build DNAscan pipeline
############################################################


FROM ubuntu:latest
MAINTAINER Heather Marriott "heather.marriott@kcl.ac.uk"
SHELL ["/bin/bash", "-c"]
RUN apt-get update
RUN apt-get install -y git
RUN cd $PWD
RUN git clone https://github.com/hevmarriott/DNAscanv2.git
RUN cd DNAscanv2
RUN bash DNAscanv2/scripts/install_dependencies_hg19.sh $PWD/DNAscanv2/local/ $PWD/DNAscanv2/ 
RUN source ~/.bashrc

WORKDIR $PWD/DNAscanv2/

ENTRYPOINT cat $PWD/DNAscanv2/docker/welcome_message.txt
