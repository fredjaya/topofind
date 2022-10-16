nextflow.enable.dsl = 2

include {

    i1 as i1_1;
    mast;
    i1 as i1_2a;
    i1 as i1_2b;
    i2

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
    mast(params.prefix, params.aln_ch, params.aln_format, i1_1.out.t2, params.nthreads)
    i1_2a(params.prefix, mast.out.class_1, "fasta", params.nthreads)
    i1_2b(params.prefix, mast.out.class_2, "fasta", params.nthreads)

    new_trees = i1_2a.out.t2
        .mix(i1_2b.out.t2)
        .flatten()

    i2(i1_2a.out.t1, i1_2b.out.t1, new_trees)
    // Count number of trees output by i1_2a + i1_2b
    // .collect(A1, A2, B1, B2).size()
    // if n_trees == 2 : TERMINATE
    // elif n_trees == 3 : mast()
    // elif n_trees == 4:
    //      .collect(A1, A2, B1).set{A1A2B1}
    //      .collect(A1, B1, B2).set{A1A2B1}

}
