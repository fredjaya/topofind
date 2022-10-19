nextflow.enable.dsl = 2

include {

    split_aln as A01_split_aln;
    mast as A02_mast_t2;
    split_aln as B01_split_aln;
    split_aln as C01__split_aln;
    sort_trees;
    mast as mast_a1b12;
    mast as mast_a12b1;

} from "./workflows.nf"

params.aln_ch = Channel
    .fromPath(params.aln)

params.aln_name = Channel
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

    // Debugging
    previous_model = ${params.previous_model}
    """
  
    A01_split_aln(params.aln_name, "A01_split_aln",  params.aln_ch, params.aln_format, params.nthreads)
   
 /*
    A02_mast_t2(params.prefix, params.aln_ch, params.aln_format, i1_1.out.t2, params.nthreads)
    B01_split_aln(params.prefix, mast_1.out.class_1, "fasta", params.nthreads)
    C01_split_aln(params.prefix, mast_1.out.class_2, "fasta", params.nthreads)
   */
 
    /*
     * Mixing channels should be conducted outside of subworkflows to ensure
     * all outputs are collected prior to downstream processes 
     */
/*
    i1_2a.out.t2
        .mix(i1_2b.out.t2)
        .flatten()
        .branch {
            class_1_split: it =~ /class_1-out_class/
            class_2_split: it =~ /class_2-out_class/
        } .set { new_trees } 

    i1_2a.out.t1.mix(new_trees.class_2_split).set { trees_a1b12 }
    i1_2b.out.t1.mix(new_trees.class_1_split).set { trees_a12b1 }

    trees_a1b12.view()
*/
    //sort_trees(i1_2a.out.t1, i1_2b.out.t1, new_trees)
    //mast_a1b12(params.prefix, params.aln_ch, params.aln_format, sort_trees.out.trees_a1b12, params.nthreads)
    //mast_a12b1(params.prefix, params.aln_ch, params.aln_format, sort_trees.out.trees_a12b1, params.nthreads)
    
    // Count number of trees output by i1_2a + i1_2b
    // .collect(A1, A2, B1, B2).size()
    // if n_trees == 2 : TERMINATE
    // elif n_trees == 3 : mast()
    // elif n_trees == 4:
    //      .collect(A1, A2, B1).set{A1A2B1}
    //      .collect(A1, B1, B2).set{A1A2B1}

    // Select the best performing three-tree mast 
    // Then, split all n blocks?

}
