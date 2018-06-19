# Installation

## Prerequisites

### [OpenMPI](http://www.open-mpi.org/)

Either install a recent OpenMPI version with your package manager or install the latest version manually. This can be done, e.g. for version 1.6.5, by:

```
wget http://www.open-mpi.org/software/ompi/v1.6/downloads/openmpi-1.6.5.tar.bz2
tar -xjvf openmpi-1.6.5.tar.bz2
cd openmpi-1.6.5
CFLAGS="-O2 -pipe -march=core2" ./configure -prefix=$HOME/chroot --disable-mpi-fortran --enable-shared --enable-static --enable-cxx-exceptions --with-sge --without-openib --enable-memchecker --enable-debug
make install
```

This can be done without root privileges and installs OpenMPI in /home/$USER/chroot

### [FFTW](http://www.fftw.org/)

Fargo needs version 3 of the FFTW library with MPI support! Most distributions only provide a non MPI version package. The latest version can be installed manually by:

```
wget http://www.fftw.org/fftw-3.3.3.tar.gz
tar -xvzf fftw-3.3.3.tar.gz   
cd fftw-3.3.3 
MPICC="/home/$(USER)/chroot/bin/mpicc" LDFLAGS="-Wl,--rpath -Wl,/home/$USER/chroot/lib" CFLAGS="-O3 -pipe -march=core2 -L/home/$USER/chroot/lib" ./configure --disable-fortran --enable-mpi  --enable-shared --prefix=$HOME/chroot
make -j
make install
```

This installs FFTW also in /home/$USER/chroot and assumes that you've installed OpenMPI there as well. If not, you need to adjust the configure line.

## FARGO

Download the latest version of the Code e.g. by
```
git clone git@github.com:twam/FARGO.git
```

Adjust the `makefile` in the `src` folder to your needs. Primarily you need to set the paths to the MPI Compilers and the Prefix Paths for MPI and FFTW. If you installed OpenMPI and FFTW as explained above you'll have to set
```
CXX = /home/$USER/chroot/bin/mpic++
MPI_PREFIX = /home/$(USER)/chroot
FFTW_PREFIX= /home/$(USER)/chroot
```

Afterwards you can compile the code by
```
make -j
```

# Use

To run a simulation create a new Folder, e.g. ~/test and create two directories in it: in and out1:

```
mkdir ~/test
cd ~/test
mkdir in
mkdir out1
```
Copy the `example.par` out of the src folder into the in folder and adjust it to your needs.
```
cp ~/FARGO/src/exampler.par in/in.par
nano in/in.par
```
Run the simulation by
```
 ~/chroot/bin/mpirun -np 4  --prefix ~/chroot ~/FARGO/fargo -v  in/in.par
```
