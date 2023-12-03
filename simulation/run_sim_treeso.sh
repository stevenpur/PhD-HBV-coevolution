#!/bin/bash
#SBATCH -D ./
#SBATCH -A bag.prj
#SBATCH -p short
#SBATCH --cpus-per-task=10
#SBATCH --array=1-100
#SBATCH --output=result-%j.out

# Define n (sample size) and f values
n_values=(100 1000 2000) # Sample values
f_values=(2 10 100)    # Coev factors values

# The run number for this job, based on the Slurm array task ID
run_number=$SLURM_ARRAY_TASK_ID

# Calculate total combinations of n and f
total_combinations=$(( ${#n_values[@]} * ${#f_values[@]} ))


# Iterate over all combinations of n and f, running the simulation once for each
for n in "${n_values[@]}"; do
    for f in "${f_values[@]}"; do
        runid=l100n${n}f${f}_${run_number}
        Rscript /users/bag/hlq763/hbv_covar3/github/simulation/treeso_simtest.r $runid
    done
done


