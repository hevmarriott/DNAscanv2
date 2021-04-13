############################################################
# Dockerfile to build DNAscan pipeline
############################################################


FROM ubuntu:latest
MAINTAINER Alfredo Iacoangeli "alfredo.iacoangeli@kcl.ac.uk"

RUN apt-get update
RUN apt-get install -y git
RUN cd ~
RUN git clone https://github.com/hevmarriott/DNAscan.git
RUN cd DNAscan
RUN bash DNAscan/scripts/install_dependencies_hg19_docker.sh ~ ~/DNAscan/ 
RUN source ~/.bashrc

WORKDIR ~/DNAscan/

ENTRYPOINT cat ~/DNAscan/docker/welcome_message.txt
