# Self-gravity solver test azimuthal

Test of the self-gravity solver in the azimuthal direction.

The density is initialized with two guassian peaks and the SG acceleration is evaluated on an azimuthal cut through the centers.

To run this test, you have to modify the source and then run `run_auto_test.sh`
At the moment, this is not implemented to run automatically.
Add
```c++
selfgravity::compute(data, 0, false);
sim::handle_outputs(data);
```
just before 
```c++
sim::run(data);	
```
in `main.cpp`.

This will recalculate the SG acceleration and output the result just after restarting from the modified initial snapshot.