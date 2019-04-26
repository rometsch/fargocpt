#!/usr/bin/env bash
# test whether outputs are the exact same after reading and writing output files
if [ "$#" -eq 1 ]; then
	N=$1
else
	N=1
fi
ODIR=outputs/binary_io_test
SDIR=test/binary_io
# change number of total output steps to selected N
./Tools/chprm.py $SDIR/fargo.par Ntot $N
rm -rf $ODIR
mpirun ./fargo $SDIR/fargo.par
cp $ODIR/gasdens$N.dat $ODIR/gasdens$N.dat.bak
mpirun ./fargo $SDIR/fargo.par -r $N
# compare the results
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "comparing results (should be empty for identical files)"
cmp -b $ODIR/gasdens$N.dat.bak $ODIR/gasdens$N.dat
echo "do statistics"
$SDIR/print_differences.py $ODIR/gasdens$N.dat.bak $ODIR/gasdens$N.dat
