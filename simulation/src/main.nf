#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define parameters
params.n_values = [400]
params.f_values = [8]
params.m_values = [1e-04, 0.001]
params.pop_sizes = [400]
params.run_indices = 1..10
params.len = 100
params.covar_desc_per = 0.5
params.outdir = "/users/bag/hlq763/hbv_covar3/analysis/simulation/trial_latest/"


// Process to run simulations
process tree_simulation {
    publishDir "${params.outdir}", mode: 'copy'
    
    // Add SLURM configuration
    executor 'slurm'
    clusterOptions '-p short'
    memory '4GB'
    cpus 1
    time '1h'
    
    input:
    tuple val(pop_size), val(run_index)

    output:
    file '*.tree'
    
    script:
    """
    #!/bin/bash
    set -e
    
    echo "Starting process for pop_size=${pop_size}, run_index=${run_index}"
    Rscript /users/bag/hlq763/hbv_covar3/github/simulation/src/sim_tree.r \
        --pop ${pop_size} \
        --runs ${run_index} \
        --rescale \
    """
    
}
// Process to run simulations
process run_simulation {
    publishDir "${params.outdir}/${param_set}", mode: 'copy'
    
    // Add SLURM configuration
    executor 'slurm'
    clusterOptions '-p short'
    memory '4GB'
    cpus 1
    time '1h'
    
    input:
    tuple val(f), val(m), val(len), val(pop_size), val(run_index), val(param_set)

    output:
    file '*.fasta'
    file '*.log'
    file tree_file

    script:
    """
    #!/bin/bash
    set -e
    
    echo "Starting process for run_index=${run_index}, f=${f}, m=${m}, len=${len}"
    log_file="simseq_${param_set}_${run_index}.log"
    
    Rscript /users/bag/hlq763/hbv_covar3/github/simulation/src/simulate_coev_seq.r \
        --len ${len} \
        --tree ${tree_file} \
        --coev_factor ${f} \
        --run_ind ${run_index} \
        --covar_desc_per ${params.covar_desc_per} \
        --mu ${m} \
        --output_dir \$PWD \
        > \${log_file} 2>&1
    """
}

// Create channels for parameters
Channel.from(params.run_indices) 
    .set { run_indices_ch }

// Create a channel with all parameter combinations
Channel.from(params.pop_sizes) 
    .combine(run_indices_ch)
    .map { pop_size, run_index -> 
        tuple(pop_size, run_index) 
    }
    .set { pop_runid_ch }

Channel.from(params.f_values)
    .combine(Channel.from(params.m_values))
    .combine(Channel.from(params.len))
    .combine(Channel.from(params.pop_sizes))
    .combine(run_indices_ch)
    .map { f, m, len, pop_size, run_index -> 
        def param_set = "l${len}n${pop_size}f${f}u${m}"
        log.info "Created parameter combination: f=${f}, m=${m}, len=${len}, pop_size=${pop_size}"
        tuple(f, m, len, pop_size, run_index, param_set) 
    }
    .set {param_combinations_ch}


// Workflow definition
workflow {
    tree_simulation(pop_runid_ch)
    run_simulation(param_combinations_ch)
} 