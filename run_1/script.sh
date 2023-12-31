#!/bin/bash
#
#PBS -N PDB_CODE_xxx
#PBS -o xxx.log
#PBS -e xxx.err
#PBS -l mem=1GB
#PBS -l walltime=6:00:00

# Calculate the number of processors allocated to this run.
NPROCS=`wc -l < $PBS_NODEFILE`

# Calculate the number of nodes allocated.
NNODES=`uniq $PBS_NODEFILE | wc -l`

### Display the job context
echo Running on host `hostname`
echo Time is `date`


cd $PBS_O_WORKDIR
icpc -o PDB_CODE_xxx HIO_boundary_hist_match.cpp -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -O -no-multibyte-chars;
./PDB_CODE_xxx
