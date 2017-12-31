#!/bin/bash
#SBATCH -J MM10MC14
#SBATCH -n 1 
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH -p zorana # Partition to submit to
#SBATCH --array=0-9
#SBATCH --mem=200 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o jobarray.out # File to which STDOUT will be written
#SBATCH -e jobarray.err # File to which STDERR will be written
#SBATCH --mail-type=END # Type of email notification- BEGIN,END,FAIL,A


./a.out $SLURM_ARRAY_TASK_ID 0425_$SLURM_ARRAY_TASK_ID
