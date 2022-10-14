nextflow.enable.dsl = 2

include {

    i1 as i1_1;
    mast;
    i1 as i1_2;

} from "./workflows.nf"

params.aln_ch = Channel
    .fromPath(params.aln)

params.prefix = Channel
    .fromPath(params.aln)
    .map { file -> file.simpleName }

workflow {
    
    log.info"""
    // Paths and directories
    base        = ${baseDir}
    out         = ${params.out}

    // Resource usage
    nthreads    = ${params.nthreads}
    ncpus       = ${params.ncpus}

    // Input alignment
    aln         = ${params.aln}
    aln_format  = ${params.aln_format}
    """
    
    i1_1(params.prefix, params.aln_ch, params.aln_format, params.nthreads)
    mast(params.prefix, params.aln_ch, params.aln_format, i1_1.out.trees, params.nthreads)
    mast.out.splitted_aln.view() 
    mast.out.bic.view() 

}
