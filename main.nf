nextflow.enable.dsl = 2

include {

    split_aln as A01_split_aln;
    mast as BC02_mast;
    split_aln as B03_split_aln;
    split_aln as C03_split_aln;
    mast as B1C1204_mast;
    mast as B12C104_mast;

} from "./workflows.nf"

include { compare_bic } from "./processes.nf"

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
    A01_split_aln.out.bic.map{ it -> it.tokenize(" ") }.set{ bic_1 }

    BC02_mast(params.aln_name, "02_BC_mast", params.aln_ch, params.aln_format, A01_split_aln.out.t2, "GTR+FO+G,GTR+FO+G", params.nthreads)
    BC02_mast.out.bic.map{ it -> it.tokenize(" ") }.set{ bic_2 }

    B03_split_aln(params.aln_name, "03_B_split_aln", BC02_mast.out.class_1, "fasta", params.nthreads)
    C03_split_aln(params.aln_name, "03_C_split_aln", BC02_mast.out.class_2, "fasta", params.nthreads)
 
    /*
     * Mixing channels should be conducted outside of subworkflows to ensure
     * all outputs are collected prior to downstream processes 
     */

    B03_split_aln.out.t1.combine(C03_split_aln.out.t2).set { trees_B1C12 }
    C03_split_aln.out.t1.combine(B03_split_aln.out.t2).set { trees_B12C1 }

    // Count number of trees output by i1_2a + i1_2b
    // .collect(A1, A2, B1, B2).size()
    // if n_trees == 2 : TERMINATE
    // elif n_trees == 3 : mast()
    // elif n_trees == 4:
    //      .collect(A1, A2, B1).set{A1A2B1}
    //      .collect(A1, B1, B2).set{A1A2B1}

    B1C1204_mast(params.aln_name, "04_B1C12_mast", params.aln_ch, params.aln_format, trees_B1C12, "GTR+FO+G,GTR+FO+G,GTR+FO+G", params.nthreads)
    B12C104_mast(params.aln_name, "04_B12C1_mast", params.aln_ch, params.aln_format, trees_B12C1, "GTR+FO+G,GTR+FO+G,GTR+FO+G", params.nthreads)

    B1C1204_mast.out.bic.map { it -> it.tokenize(" ") }.set{ bic_3a }
    B12C104_mast.out.bic.map { it -> it.tokenize(" ") }.set{ bic_3b }
  
    // Combine BIC scores in a single data structure?
    bic_1
        .concat(bic_2)
        .concat(bic_3a)
        .concat(bic_3b)
        .collect()
        .map { it -> it.collate(2) } 
        .set { bic_all }

    bic_all.view()

    // Select the best performing three-tree mast 
    
    x = bic_all.map { it -> it[2] }
    y = bic_all.map { it -> it[3] }

    compare_bic(x, y)
}
