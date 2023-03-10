nextflow.enable.dsl = 2

include {
    split_aln;
    mast;
} from "./workflows.nf"

include {
    recurse_trees;
} from "./processes.nf"

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

    n_trees = 3
    run_names = ["3_B_AA_AB", "3_A_BA_BB"]
    if (n_trees == 1) {
        run_names = '2_A_B'
    } else {
        println recurse_trees(run_names) 
    }
    
    //split_aln("01_split_A_B", params.aln_ch, params.nthreads)
    //mast("02_mast_A_B", params.aln_ch, split_aln.out.t2, "GTR+FO+G,GTR+FO+G", params.nthreads)
}
