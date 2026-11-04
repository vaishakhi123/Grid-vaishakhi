# OLCF Riker - GPU partition: 8 nodes, each 1x AMD EPYC 9575F (64 cores, 2 NUMA domains)
# + 2x NVIDIA L40S (Ada, sm_89, 48 GB, PCIe, no NVLink), InfiniBand.
# https://docs.olcf.ornl.gov/systems/riker_user_guide.html
# GPU-aware MPI (per the guide): cuda/13.3.0 must be loaded before mpich/5.0.1-gpu,
# and MPIR_CVAR_ENABLE_GPU=1 at run time.
module reset >/dev/null 2>&1
module load gcc/14.2.0
module load cuda/13.3.0
module load mpich/5.0.1-gpu
module load cmake/3.31.11 gmp/6.3.0 fftw/3.3.11-omp hdf5/1.14.6

export MPIR_CVAR_ENABLE_GPU=1

# install prefix for Grid, Hadrons and the small dependencies built from source
export RIKER_SW=${RIKER_SW:-$HOME/riker-sw}
export DEPS=$RIKER_SW/deps

# module install prefixes: spack modules prepend the prefix to CMAKE_PREFIX_PATH
module_root() { module show "$1" 2>&1 | grep CMAKE_PREFIX_PATH | grep -oE '/[^",)]+' | head -1; }
export CUDA_HOME=${CUDA_HOME:-$(dirname $(dirname $(which nvcc)))}
export GMP=${GMP:-$(module_root gmp)}
export FFTW=${FFTW:-$(module_root fftw)}
export HDF5=${HDF5:-$(module_root hdf5)}
export CLIME=$DEPS
export MPFR=$DEPS
# hdf5 module without the C++ API -> build_deps.sh puts one in $DEPS
[ -f $HDF5/include/H5Cpp.h ] || export HDF5=$DEPS

export LD_LIBRARY_PATH=$DEPS/lib:$GMP/lib:$FFTW/lib:$HDF5/lib:$RIKER_SW/grid/lib:$LD_LIBRARY_PATH
unset CC CXX MPICXX MPICC
ulimit -c 0
