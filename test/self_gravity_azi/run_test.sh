#!/usr/bin/env bash
cd $(dirname $0)
fargocpt run start -N 0 setup.yml # > /dev/null
python3 custom_init.py ../../output/tests/self_gravity_solver_azi/out
# rm ../../output/tests/self_gravity_solver_2D/out/snapshots/0/a_sg*
fargocpt run auto setup.yml # > /dev/null
python3 ../run_test.py --silent -t