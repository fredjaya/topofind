nextflow.enable.dsl = 2

include {

    split_aln as A01_split_aln;
    mast as BC02_mast;
    split_aln as B03_split_aln;
    split_aln as C03_split_aln;
    mast as B1C1204_mast;
    mast as B12C104_mast;

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
  
    A01_split_aln(params.aln_name, "01_A_split_aln",  params.aln_ch, params.aln_format, params.nthreads)
    BC02_mast(params.aln_name, "02_BC_mast", params.aln_ch, params.aln_format, A01_split_aln.out.t2, params.nthreads)
    B03_split_aln(params.aln_name, "03_B_split_aln", BC02_mast.out.class_1, "fasta", params.nthreads)
    C03_split_aln(params.aln_name, "03_C_split_aln", BC02_mast.out.class_2, "fasta", params.nthreads)
 
    /*
     * Mixing channels should be conducted outside of subworkflows to ensure
     * all outputs are collected prior to downstream processes 
     */
    B03_split_aln.out.t2
        .mix(C03_split_aln.out.t2)
        .flatten()
        .branch {
            block_B_split: it =~ /class_1-out_class/
            block_C_split: it =~ /class_2-out_class/
        } .set { new_trees } 

    B03_split_aln.out.t1.mix(new_trees.block_B_split).set { trees_B1C12 }
    C03_split_aln.out.t1.mix(new_trees.block_C_split).set { trees_B12C1 }

    trees_B1C12.view()

    B1C1204_mast(params.aln_name, "04_B1C12_mast", params.aln_ch, params.aln_format, trees_B1C12, params.nthreads)
    B12C104_mast(params.aln_name, "04_B12C1_mast", params.aln_ch, params.aln_format, trees_B12C1, params.nthreads)
    
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
