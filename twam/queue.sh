#!/usr/bin/env bash

source sheep_env

case $(hostname -s) in
	login0*)
		qsub -N "$(cat uuid)" \
			 -e err.log \
			 -o out.log \
			 -m bae \
			 -M "thomas.rometsch@uni-tuebingen.de" \
			 -l walltime=03:00:00 \
			 -l nodes=1:ppn=$NCPU \
			 -q $QUEUE \
			 $1
		;;
	*)
	    $PWD/$1 1>out.log 2>err.log &
		;;
esac
