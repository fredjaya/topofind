nextflow.enable.dsl = 2
nextflow.preview.recursion=true

include {
    test;
    split_aln;
    mast;
} from "./workflows.nf"

include { update_run_names } from "./processes.nf"

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
    
    //run_names = "['3_A_BA_BB', '3_B_AA_AB']"
    //run_names = "null"
    //test.recurse(run_names, file(params.aln_ch), params.nthreads).times(2)
    //test(run_names, params.aln_ch, params.nthreads)
    
    run_name = "null"
    update_run_names(run_name) 
    new_runs = update_run_names.out.map{ it -> it.trim() }
    split_aln(new_runs, params.aln_ch, params.nthreads)
    mast(new_runs, params.aln_ch, split_aln.out.PartitionedTrees, "GTR+FO+G,GTR+FO+G", params.nthreads)
    // TODO: Add BIC comparison here
    //mast.out.aln.view()
}
