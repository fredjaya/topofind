nextflow.enable.dsl = 2
nextflow.preview.recursion=true

include { iterative } from "./workflows.nf"
include { update_run_names } from "./processes.nf"

params.aln_ch = Channel
    .fromPath(params.aln)

params.aln_name = Channel
    .fromPath(params.aln)
    .map { file -> file.simpleName }

workflow {

    log.info"""
    base        = ${baseDir}
    out         = ${params.out}
    aln         = ${params.aln}
    nthreads    = ${params.nthreads}
    """
    
    run_name = "null"
    iterative(run_name, params.aln_ch, params.nthreads)
}
