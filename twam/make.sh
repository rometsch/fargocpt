#!/usr/bin/env bash

case $(hostname -s) in
	login0*)
		module add mpi/openmpi/1.10.3-gnu-5.2
		module add numlib/fftw/3.3.5-openmpi-1.10.3-gnu-5.2
		module add numlib/gsl/2.1
		export FARGO_ARCH="BINAC"
		;;
	cpt-titan)
		export FARGO_ARCH="CPTTITAN"
		;;
	cpt-pandora)
		export FARGO_ARCH="LINUX"
		;;
	cpt-*)
		export FARGO_ARCH="CPT"
		;;
  	*)
		export FARGO_ARCH="LINUX"
  		;;
esac

cd src
make
cd ..
