#!/usr/bin/env bash

case $(hostname -s) in
	login0*)
		# We are on binac
		module load devel/cuda/8.0
		qsub -N "$(cat uuid)" \
			 -e err.log \
			 -o out.log \
			 -m bae \
			 -M "thomas.rometsch@uni-tuebingen.de" \
			 -l walltime=03:00:00 \
			 -l nodes=1:ppn=1:gpus=1:exclusive_process \
			 -q gpu \
			 $1
		;;
	cpt-titan)
		bsub -o out.log -e err.log -q default $PWD/$1
		;;

	cpt-*)
		bsub -o out.log -e err.log -q gpu0 $PWD/$1
		;;
  	*)
  		;;
esac
