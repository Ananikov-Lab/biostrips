FROM python:3.9

RUN mkdir /src
WORKDIR /src
COPY . /src
RUN apt-get update
RUN pip install -r requirements.txt

