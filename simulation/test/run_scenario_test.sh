#!/bin/sh

# run the simulaation at different mutation rates
for mu in 0 0.01 0.02 0.04 0.08
    do
        # run the simulation with the given mutation rate
        echo "Running simulation with mutation rate: $mu"
        
         Rscript $HOME/hbv_covar3/github/simulation/simulate_coev_seq_scenario_test.r \
             --len 100 \
             --tree $HOME/hbv_covar3/analysis/sim_seq/simseq_N1600_1_rescaled.tree \
             --coev_factor 16 \
             --run_ind 1 \
             --mu ${mu} \
             --covar_desc_per 0.5 \
             --out $HOME/hbv_covar3/analysis/sim_seq/test/test_comut_n1600u${mu}f16per0.5.fasta
        # run the pydca analysis on it
        echo "Running pydca analysis on simulated sequence"
        python $HOME/hbv_covar3/github/simulation/pydca_test_scenario.py \
            $HOME/hbv_covar3/analysis/sim_seq/test/test_comut_n1600u${mu}f16per0.5.fasta
    done


