FROM continuumio/miniconda3

RUN mkdir /src
WORKDIR /src
COPY . /src
RUN apt-get update
RUN conda env create -f env.yml
RUN conda init bash
RUN echo "conda activate biostrips" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]

