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
f_values=(2 4 8 16 32 64 128 )    # Coev factors values
m_values=(1e-04 0.001 0.01 0.05 0.25 1.25 6.25)

# The run number for this job, based on the Slurm array task ID
# run_number=$SLURM_ARRAY_TASK_ID

# hard code run number to 1
run_number=1
# Calculate total combinations of n and f
total_combinations=$(( ${#n_values[@]} * ${#f_values[@]} ))

# get into the right directory
cd /well/bag/clme1992/hbv_covar3/analysis/sim_seq

# Iterate over all combinations of n and f, running the simulation once for each
for n in "${n_values[@]}"; do
    tree_file="../simseq_N${n}_${run_number}_rescaled.tree"
    for f in "${f_values[@]}"; do
        for m in "${m_values[@]}"; do
            param_set="l100n${n}f${f}u${m}comut0.5"
            echo "simulating sequence data for ${param_set}..."
            mkdir -p ${param_set}_test
            cd ${param_set}_test
            log_file="simseq_${param_set}_${run_number}_rescaled.log"
            Rscript /users/bag/hlq763/hbv_covar3/github/simulation/simulate_coev_seq_scenario_test.r \
                --len 100 \
                --tree ${tree_file} \
                --coev_factor ${f} \
                --run_ind $run_number \
                --covar_desc_per 0.5 \
                --mu ${m} > ${log_file} 2>&1 \
                --out simseq_${param_set}_${run_number}_rescaled.fasta
            cd ..
            echo "simulation complete!"
        done
    done
done


