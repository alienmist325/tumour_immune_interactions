#!/bin/bash
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -l walltime=00:01:00
#PBS -N numpy
 
cd $PBS_O_WORKDIR
 
module load tools/prod
module load Python/3.11.2-GCCcore-12.2.0-bare
source sim_venv/bin/activate

python tests/numpy_test.py > log.txt