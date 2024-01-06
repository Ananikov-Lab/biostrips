FROM continuumio/miniconda3

RUN mkdir /src
WORKDIR /src
COPY . /src
RUN apt-get update
RUN conda env create -f env.yml
RUN conda init bash
SHELL ["conda", "run", "-n", "biostrips", "/bin/bash", "-c"]

