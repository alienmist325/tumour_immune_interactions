#!/bin/bash
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -l walltime=00:01:00
#PBS -N hello_world
 
cd $PBS_O_WORKDIR
 
module load tools/prod
module load SciPy-bundle/2022.05-foss-2022a

# cp hello_world.py $TMPDIR
# cd $TMPDIR

python hello_world.py > log.txt
cp * $HOME/tumour_immune_interactions/job_data