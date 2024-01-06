FROM continuumio/miniconda3

RUN mkdir /src
WORKDIR /src
COPY . /src
RUN apt-get update
RUN conda env create -f env.yml
RUN conda init bash
RUN echo "source activate biostrips" > ~/.bashrc
ENV PATH /opt/conda/envs/biostrips/bin:$PATH
