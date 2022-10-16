nextflow.enable.dsl = 2

include {

    i1 as i1_1;
    mast as mast_1;
    i1 as i1_2a;
    i1 as i1_2b;
    sort_trees;
    mast as mast_a1b12;
    mast as mast_a12b1;


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
    mast_1(params.prefix, params.aln_ch, params.aln_format, i1_1.out.t2, params.nthreads)
    i1_2a(params.prefix, mast_1.out.class_1, "fasta", params.nthreads)
    i1_2b(params.prefix, mast_1.out.class_2, "fasta", params.nthreads)

    new_trees = i1_2a.out.t2
        .mix(i1_2b.out.t2)
        .flatten()

    sort_trees(i1_2a.out.t1, i1_2b.out.t1, new_trees)
    mast_a1b12(params.prefix, params.aln_ch, params.aln_format, sort_trees.out.trees_a1b12, params.nthreads)
    mast_a12b1(params.prefix, params.aln_ch, params.aln_format, sort_trees.out.trees_a12b1, params.nthreads)
    // Count number of trees output by i1_2a + i1_2b
    // .collect(A1, A2, B1, B2).size()
    // if n_trees == 2 : TERMINATE
    // elif n_trees == 3 : mast()
    // elif n_trees == 4:
    //      .collect(A1, A2, B1).set{A1A2B1}
    //      .collect(A1, B1, B2).set{A1A2B1}

}
