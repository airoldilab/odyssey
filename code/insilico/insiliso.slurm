#! /bin/sh
#SBATCH -n 1               # (Max) number of tasks per job, for R usually 1
#SBATCH -o out-%a.txt      # File for the standard output
#SBATCH -e err-%a.txt      # File for the standard error
#SBATCH -p general         # Partition to use
#SBATCH --mem-per-cpu=1024 # Memory required per CPU, in MegaBytes
#SBATCH -a 1-500           # Array of 500 jobs, with ids 1, 2, ..., 500

# And the code to run:
Rscript -e 'library(R221) ; nnmRunInsilico()' ${SLURM_ARRAY_TASK_ID}
