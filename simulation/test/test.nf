#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define parameters
params.n_values = [400, 800, 1600]
params.f_values = [8, 16, 32]

// Create channels for parameters
Channel.from(params.n_values) 
    .combine(Channel.from(params.f_values))
    .map { n, f -> 
        tuple(n, f)
    }
    .set { n_values_ch }

// Process to run simulations
process run_simulation {
    publishDir "/users/bag/hlq763/hbv_covar3/analysis/simulation/test", mode: 'copy'
    
    // Add SLURM configuration
    // executor 'slurm'
    // clusterOptions '-p short'
    // memory '4GB'
    // cpus 1
    // time '1h'
    
    input:
    tuple val(n), val(f) 

    output:
    file ("${n}_${f}.log")
    
    script:
    """
    echo "Starting process for n=${n}, f=${f}" > ${n}_${f}.log
    """
    
}

// Workflow definition
workflow {
    run_simulation(n_values_ch)
} 