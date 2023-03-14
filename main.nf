nextflow.enable.dsl = 2

include { iteration } from "./workflows.nf"

params.base_aln_ch = Channel.fromPath(params.base_aln)
params.part_aln_ch = Channel.fromPath(params.part_aln)

workflow {

    /*
     * Set default parameters
     */

    log.info"""
    base_dir    = ${baseDir}
    out_dir     = ${params.out_dir}
    base_aln    = ${params.base_aln}
    part_aln    = ${params.part_aln}
    nthreads    = ${params.nthreads}
    prev_runs   = ${params.prev_runs}
    """
   
    /*
     * Empty variables for the first iteration
     */ 
    iteration(params.prev_runs, params.base_aln_ch, params.part_aln_ch, params.nthreads)

}
