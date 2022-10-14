nextflow.enable.dsl = 2

params.aln_ch = Channel
    .fromPath(params.aln)

params.prefix = Channel
    .fromPath(params.aln)
    .map { file -> file.simpleName }

def store_models(x) {
    // Store number of FreeRate categories, model, and BIC scores
    x.tokenize(" ").collate(2)
}

def compare_bic(models_list, rk1, rk2) {
    bic1 = models_list[rk1][2]
    bic2 = models_list[rk2][2]

    if( bic2 < bic1 )
        return true
    else
        return false
}

include { 
    t1_modelfinder_across_rhas_categories; 
    keep_modelfinder_results; best_model_per_rhas; 
    t1_iqtree_with_best_r2;
    get_bic as get_bic_t1_r2;
    hmm_assign_sites as hmm_assign_sites_t1_r2; 
    evaluate_partitions; 
    split_aln;
    t1_iqtree_per_split;
    get_bic as get_bic_t1_split;
    concatenate_trees_for_mast; 
    t2_iqtree_mast; 
    get_bic as get_bic_t2_mast;
    hmm_assign_sites as hmm_assign_sites_mast;
} from './processes.nf' 

workflow i1 {
    
    /*
     * Split an alignment into two blocks according to R2 FreeRate classes
     *
     */
    take:
        prefix
        aln_ch
        aln_format
        nthreads

    main:
        t1_modelfinder_across_rhas_categories(prefix, aln_ch, nthreads)
        keep_modelfinder_results(prefix, t1_modelfinder_across_rhas_categories.out[2])
        best_model_per_rhas(prefix, keep_modelfinder_results.out)
        best_model_per_rhas.out[2]
            .map { models_out -> store_models(models_out) }
            .set { models_list }
        t1_iqtree_with_best_r2(prefix, aln_ch, nthreads, models_list.map { x -> x[1][0] } )
        get_bic_t1_r2(prefix, t1_iqtree_with_best_r2.out[2], "t1_r2")
        hmm_assign_sites_t1_r2(prefix, t1_iqtree_with_best_r2.out[5], t1_iqtree_with_best_r2.out[4], "t1_r2")
        evaluate_partitions(prefix, hmm_assign_sites_t1_r2.out[1])
        split_aln(prefix, aln_ch, evaluate_partitions.out[0], aln_format)
        t1_iqtree_per_split(split_aln.out[0].flatten(), nthreads)

    emit:
        trees = t1_iqtree_per_split.out[1].collect()
}

workflow mast {
    
    /*
     * Collect all input trees, run MAST +TR, record the BIC and 
     * assign sites with the HMM
     */

    take:
        prefix
        aln_ch
        trees
        nthreads

    main:
        concatenate_trees_for_mast(prefix, aln_ch, trees, "class_1_2")
        t2_iqtree_mast(prefix, aln_ch, nthreads, concatenate_trees_for_mast.out[0])
        get_bic_t2_mast(prefix, t2_iqtree_mast.out[2], "mast_tr_r2")
        hmm_assign_sites_mast(prefix, t2_iqtree_mast.out[4], t2_iqtree_mast.out[6], "mast_tr")

}

workflow i2 {

    /*
     * Something like run on three trees
     */

    evaluate_partitions(params.prefix, hmm_assign_sites_t1_r2.out[1])
    split_aln(params.prefix, params.aln_ch, evaluate_partitions.out[0], params.aln_format)
    t1_iqtree_per_split(split_aln.out[0].flatten(), params.nthreads)
    t2_iqtree_mast(params.prefix, params.aln_ch, params.nthreads, concatenate_trees_for_mast.out[0])
    concatenate_trees_for_mast(params.prefix, params.aln_ch, t1_iqtree_per_split.out[1].collect(), "class_1_2")
    get_bic_t2_mast(params.prefix, t2_iqtree_mast.out[2], "mast_tr_r2")
    hmm_assign_sites_mast(params.prefix, t2_iqtree_mast.out[4], t2_iqtree_mast.out[6], "mast_tr")

}

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
    
    i1(params.prefix, params.aln_ch, params.aln_format, params.nthreads)
    mast(params.prefix, params.aln_ch, i1.out.trees, params.nthreads)

}
