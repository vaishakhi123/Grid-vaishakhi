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

source $HOME/riker/Riker/sourceme.sh
SEL=$HOME/riker/Riker/select_gpu
BENCH=$RIKER_SW/src/Grid/build/benchmarks/Benchmark_staggered
nvidia-smi topo -m

srun -n1 -c32 --gpu-bind=none $SEL $BENCH --grid 16.16.16.32 --mpi 1.1.1.1 --device-mem 16000 --shm 2048
srun -n2 -c32 --gpu-bind=none $SEL $BENCH --grid 16.16.16.32 --mpi 1.1.1.2 --device-mem 16000 --shm 2048
# if the 2-rank run fails in intra-node IPC (PCIe, no NVLink), retry with --shm-mpi 1
