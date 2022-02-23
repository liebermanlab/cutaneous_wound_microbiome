#!/bin/bash
#SBATCH -p sched_mem1TB,defq
#SBATCH -n 8
#SBATCH --time=30:00:00
#SBATCH -o masterout.txt
#SBATCH -e mastererr.txt
#SBATCH --mem=128000


bash unifrac_run.sh

echo Done!!!
