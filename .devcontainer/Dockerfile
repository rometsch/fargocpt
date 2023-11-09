# syntax=docker/dockerfile:1
FROM mcr.microsoft.com/devcontainers/universal:2

RUN apt-get update -qq
RUN apt-get install -qq -y --no-install-recommends apt-utils
RUN echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections
RUN apt-get install -qq -y -o=Dpkg::Use-Pty=0 build-essential make
RUN apt-get install -qq -y -o=Dpkg::Use-Pty=0 git
RUN apt-get install -qq -y -o=Dpkg::Use-Pty=0 libopenmpi-dev 
RUN apt-get install -qq -y -o=Dpkg::Use-Pty=0 libgsl-dev 
RUN apt-get install -qq -y -o=Dpkg::Use-Pty=0 libfftw3-mpi-dev libfftw3-dev
RUN apt-get install -qq -y -o=Dpkg::Use-Pty=0 python3 python3-pip

RUN pip install --upgrade pip
RUN pip install psutil
