# syntax=docker/dockerfile:1
FROM ubuntu:22.04

WORKDIR /app

RUN apt-get update -qq
RUN apt-get install -qq -y --no-install-recommends apt-utils
RUN echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections
RUN apt-get install -qq -y -o=Dpkg::Use-Pty=0 build-essential make
RUN apt-get install -qq -y -o=Dpkg::Use-Pty=0 git
RUN apt-get install -qq -y -o=Dpkg::Use-Pty=0 openmpi-bin libopenmpi-dev 
RUN apt-get install -qq -y -o=Dpkg::Use-Pty=0 libgsl-dev 
RUN apt-get install -qq -y -o=Dpkg::Use-Pty=0 libfftw3-mpi-dev libfftw3-dev libfftw3-3
RUN apt-get install -qq -y -o=Dpkg::Use-Pty=0 python3
RUN apt-get install -qq -y -o=Dpkg::Use-Pty=0 gosu

RUN mkdir /output
COPY .git .git
COPY src src
COPY bin/fargocpt.py bin/
COPY run_fargo .
COPY testconfig.yml .
RUN ln -s /output/ ./output

# RUN groupadd -g 999 fargouser && \
#     useradd -r -u 999 -g fargouser fargouser

# RUN chown -R fargouser:fargouser /app
# RUN chown -R fargouser:fargouser /output

# USER fargouser

RUN git reset
RUN make -C /app/src -j 4

# USER root
RUN cd /usr/local/bin && ln -s /app/run_fargo fargocpt

COPY docker/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]

CMD ["fargocpt", "-np", "1", "start", "setup.yml"]