#!/bin/bash
# Smoke test of the Grid build: staggered Dslash on 1 and 2 L40S of one node.
#SBATCH -A <project_id>
#SBATCH -J grid-stag
#SBATCH -o %x-%j.out
#SBATCH -p gpu
#SBATCH -t 00:20:00
#SBATCH -N 1
#SBATCH --gpus-per-node=2
#SBATCH --exclusive

# SRC = where build.sh put the sources (Lustre scratch); export it before sbatch or edit here
GRID_SRC=${SRC:-$HOME/riker-sw/src}/Grid
source $GRID_SRC/systems/Riker/sourceme.sh
SEL=$GRID_SRC/systems/Riker/select_gpu
BENCH=$GRID_SRC/build/benchmarks/Benchmark_staggered
nvidia-smi topo -m

srun -n1 -c32 --gpu-bind=none $SEL $BENCH --grid 16.16.16.32 --mpi 1.1.1.1 --device-mem 16000 --shm 2048
srun -n2 -c32 --gpu-bind=none $SEL $BENCH --grid 16.16.16.32 --mpi 1.1.1.2 --device-mem 16000 --shm 2048
# if the 2-rank run fails in intra-node IPC (PCIe, no NVLink), retry with --shm-mpi 1
