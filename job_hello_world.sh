#!/bin/bash
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -l walltime=00:01:00
#PBS -N hello_world
 
 
module load tools/prod

cd $PBS_O_WORKDIR
cp tests/hello_world.py $TMPDIR

cd $TMPDIR

python tests/hello_world.py > log.txt

mkdir $HOME/tumour_immune_interactions/job_data
cp log.txt $HOME/tumour_immune_interactions/job_data
cp hello_world.* $HOME/tumour_immune_interactions/job_data

cd $HOME/tumour_immune_interactions/job_data