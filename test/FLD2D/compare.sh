#!/usr/bin/env bash

FILEDIR="$(dirname $(realpath $0))"
cd $FILEDIR


REF1_NP=1
REF1_NT=4

REF2_NP=4
REF2_NT=1

rm -rf output/ref1 output/ref2

sed -i 's/RadiativeDiffusionTestModule:.*$/RadiativeDiffusionTestModule: yes/' setup.yml
sed -i 's/RadiativeDiffusionSolver:.*$/RadiativeDiffusionSolver: SOR/' setup.yml

../../run_fargo -nt 1 -N 0 start setup.yml 1> out.log 2>err.log
./create_input.py
# ../../run_fargo -nt 1 -v restart 0 setup.yml # 1> out.log 2>err.log
mpirun -x OMP_NUM_THREADS=$REF1_NT -n $REF1_NP ../../bin/fargocpt -v restart 0 setup.yml # 1> out.log 2>err.log

mv output/out output/ref1
mv out.log output/ref1
mv err.log output/ref1

sed -i 's/RadiativeDiffusionTestModule:.*$/RadiativeDiffusionTestModule: yes/' setup.yml
sed -i 's/RadiativeDiffusionSolver:.*$/RadiativeDiffusionSolver: SOR/' setup.yml

../../run_fargo -nt 1 -N 0 start setup.yml 1> out.log 2>err.log
./create_input.py
# ../../run_fargo -nt 1 -v restart 0 setup.yml # 1> out.log 2>err.log
mpirun -x OMP_NUM_THREADS=$REF2_NT -n $REF2_NP ../../bin/fargocpt -v restart 0 setup.yml 1> out.log 2>err.log

mv output/out output/ref2
mv out.log output/ref2/
mv err.log output/ref2/

diff output/ref1/snapshots/1/Trad.dat output/ref2/snapshots/1/Trad.dat

./plot_difference.py 1 output/ref1 output/ref2 -n Trad