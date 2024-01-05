#!/bin/bash
#SBATCH -D ./
#SBATCH -A bag.prj
#SBATCH -p short
#SBATCH --cpus-per-task=1
#SBATCH --array=1-5
#SBATCH --output=result-%j.out

# Define n (sample size) and f values
n_values=(100 1000 2000) # Sample values
f_values=(2 10 100)    # Coev factors values
m_values=(0.1 1 100 1000)

# The run number for this job, based on the Slurm array task ID
run_number=$SLURM_ARRAY_TASK_ID

# Calculate total combinations of n and f
total_combinations=$(( ${#n_values[@]} * ${#f_values[@]} ))

# get into the right directory
cd /well/bag/clme1992/hbv_covar3/analysis/sim_seq

# Iterate over all combinations of n and f, running the simulation once for each
for n in "${n_values[@]}"; do
    Rscript /users/bag/hlq763/hbv_covar3/github/simulation/sim_tree.r $n
    tree_file="simseq_N${n}.tree"
    for f in "${f_values[@]}"; do
        for m in "${m_values[@]}"; do
            Rscript /users/bag/hlq763/hbv_covar3/github/simulation/simulate_coev_seq.r --len 100 --tree ${tree_file} --coev_factor ${f} --run_ind $run_number --mu ${m}
        done
    done
done


