#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.iterations = 100
params.scriptDir = "/well/bag/clme1992/hbv_covar3/github/analysis"

// Create a channel emitting the numbers 1 to 100
iterations = channel.create(1..params.iterations)

process simulateCoevSeq {
    input:
    val index from iterations

    """
    Rscript ${params.scriptDir}/simulate_coev_seq.r ${index}
    """
}

process treesoSimTest {
    input:
    val index from iterations

    """
    Rscript ${params.scriptDir}/treeso_simtest.r ${index}
    """
}

workflow{
    iterations
    simulateCoevSeq(iterations)
    treesoSimTest(iterations)
}