#!/bin/bash
#SBATCH -D ./
#SBATCH -A bag.prj
#SBATCH -p short
#SBATCH --cpus-per-task=3
#SBATCH --array=1-100

source $HOME/.bashrc
source activate hbv_covar

n_values=(100 1000 2000) # Sample values
f_values=(2 10 100)    # Coev factors values
m_values=(0.01 0.1 1 10 100)

# The run number for this job, based on the Slurm array task ID
run_number=$SLURM_ARRAY_TASK_ID

# Calculate total combinations of n and f
total_combinations=$(( ${#n_values[@]} * ${#f_values[@]} * ${#m_values[@]} ))
cd /well/bag/clme1992/hbv_covar3/analysis/sim_seq/

# Iterate over all combinations of n and f, running the simulation once for each
for n in "${n_values[@]}"; do
    for f in "${f_values[@]}"; do
        for m in "${m_values[@]}"; do
            echo "running simulation for n=$n, f=$f, m=$m"
            param_id=l100n${n}f${f}u${m}
            runid=${param_id}_${run_number}
            mkdir -p $param_id
            cd $param_id 
            # log the output
            log_file=cout_${runid}.txt
            Rscript /users/bag/hlq763/hbv_covar3/github/simulation/treeso_simtest.r $runid 3 > $log_file 2>&1
            cd ..
        done
    done
done


