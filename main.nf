nextflow.enable.dsl = 2

include {

    i1 as i1_1;
    mast;
    i1 as i1_2a;
    i1 as i1_2b;

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
    i1_2a(params.prefix, mast.out.class_1, "fasta", params.nthreads)
    i1_2b(params.prefix, mast.out.class_2, "fasta", params.nthreads)

    n_trees = i1_2a.out.trees
        .mix(i1_2b.out.trees)
        .flatten()
        .count()

    // Count number of trees output by i1_2a + i1_2b
    // .collect(A1, A2, B1, B2).size()
    // if n_trees == 2 : TERMINATE
    // elif n_trees == 3 : mast()
    // elif n_trees == 4:
    //      .collect(A1, A2, B1).set{A1A2B1}
    //      .collect(A1, B1, B2).set{A1A2B1}

}
