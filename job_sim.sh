#!/bin/bash
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -l walltime=00:01:00
#PBS -N numpy
 
cd $PBS_O_WORKDIR
mkdir $TMPDIR/sim
mkdir $TMPDIR/sim_data
cp sim/* $TMPDIR/sim
cp sim_data/sim.pickle $TMPDIR/sim_data/sim.pickle


 
module load tools/prod
module load Python/3.11.2-GCCcore-12.2.0-bare
source sim_venv/bin/activate

cd $TMPDIR

python sim/main.py -sf y -ow y -c Test > log.txt

mkdir $HOME/tumour_immune_interactions/job_data
cp * $HOME/tumour_immune_interactions/job_data