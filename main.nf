nextflow.enable.dsl = 2

include { iteration } from "./workflows.nf"

workflow {

    /*
     * Set default parameters
     */

    log.info"""
    base_dir    = ${baseDir}
    out_dir     = ${params.out_dir}
    base_aln    = ${params.base_aln}
    run_file    = ${params.run_file}
    nthreads    = ${params.nthreads}
    """

    iteration(params.base_aln, params.run_file, params.submodel, params.nthreads)

}
