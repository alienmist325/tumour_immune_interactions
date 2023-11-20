#!/bin/bash
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -l walltime=00:01:00
#PBS -N numpy
 
cd $PBS_O_WORKDIR
 
module load tools/prod
module load SciPy-bundle/2022.05-foss-2022a
source sim_venv/bin/activate

python tests/numpy.py > log.txt