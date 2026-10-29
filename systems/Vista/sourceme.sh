# TACC Vista - NVIDIA Grace-Hopper GH200 (ARM + H100, sm_90)
# nvidia/24.7 (CUDA 12.5) and nvpl/24.7 are loaded by default.
# Switch to openmpi built with nvc++ for GPU-aware MPI.
module unload openmpi
module load openmpi/5.0.5_nvc249

source $HOME/spack/share/spack/setup-env.sh
spack load c-lime%gcc
spack load gmp%gcc
spack load mpfr%gcc
spack load openssl%gcc
spack load hdf5+cxx~mpi%gcc
spack load fftw%gcc

export CLIME=`spack find --paths c-lime%gcc   | grep ^c-lime | awk '{print $2}'`
export GMP=`spack find   --paths gmp%gcc      | grep ^gmp   | awk '{print $2}'`
export MPFR=`spack find  --paths mpfr%gcc     | grep ^mpfr  | awk '{print $2}'`
export HDF5=`spack find  --paths hdf5+cxx~mpi%gcc | grep ^hdf5  | awk '{print $2}'`
export FFTW=`spack find  --paths fftw%gcc     | grep ^fftw  | awk '{print $2}'`


export LD_LIBRARY_PATH=$CLIME/lib:$GMP/lib:$MPFR/lib:$HDF5/lib:$FFTW/lib:$LD_LIBRARY_PATH

unset CC CXX MPICXX MPICC

ulimit -c 0

