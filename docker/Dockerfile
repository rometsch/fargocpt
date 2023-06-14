# syntax=docker/dockerfile:1
FROM ubuntu:22.04

WORKDIR /app

RUN apt-get update -qq
RUN apt-get install -qq -y --no-install-recommends apt-utils
RUN echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections
RUN apt-get install -qq -y -o=Dpkg::Use-Pty=0 build-essential make
RUN apt-get install -qq -y -o=Dpkg::Use-Pty=0 git
RUN apt-get install -qq -y -o=Dpkg::Use-Pty=0 libopenmpi-dev 
RUN apt-get install -qq -y -o=Dpkg::Use-Pty=0 libgsl-dev 
RUN apt-get install -qq -y -o=Dpkg::Use-Pty=0 libfftw3-mpi-dev libfftw3-dev
RUN apt-get install -qq -y -o=Dpkg::Use-Pty=0 python3

COPY .git .git
COPY src src

RUN make -C /app/src -j 4

COPY bin/fargocpt.py bin/
COPY run_fargo .

RUN cd /usr/local/bin && ln -s /app/run_fargo fargocpt

WORKDIR /project

CMD ["bash"]