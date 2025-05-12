#!/bin/bash
#SBATCH -D ./
#SBATCH -A bag.prj
#SBATCH -p short
#SBATCH --cpus-per-task=1
#SBATCH --array=1-100

source $HOME/.bashrc
source activate hbv_covar

# Define n (sample size) and f values
n_values=(100 200 400 800 1600) # Sample values
m_values=(0.01 0.05 0.25 1.25 6.25 31.25 156.25)

# The run number for this job, based on the Slurm array task ID
run_number=$SLURM_ARRAY_TASK_ID

# get into the right directory
cd /well/bag/clme1992/hbv_covar3/analysis/sim_seq

# Iterate over all combinations of n and f, running the simulation once for each
for n in "${n_values[@]}"; do
    tree_file="../simseq_N${n}_${run_number}_rescaled.tree"
        for m in "${m_values[@]}"; do
            echo "simulating sequence data for n=${n},mu=${m}..."
            param_set="l100n${n}u${m}"
            mkdir -p ${param_set}
            cd ${param_set}
            log_file="simseq_${param_set}_${run_number}_rescaled.log"
            Rscript /users/bag/hlq763/hbv_covar3/github/simulation/simulate_coev_seq_null.r \
                --len 100 \
                --tree ${tree_file} \
                --run_ind $run_number \
                --mu ${m} > ${log_file} 2>&1
            cd ..
            echo "simulation complete!"
        done
done


