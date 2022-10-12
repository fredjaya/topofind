nextflow.enable.dsl = 2

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

workflow {
    
    /*
     * This workflow will only use R2 rate categories to to delimit sites 
     */ 
    
    params.aln_ch = Channel
        .fromPath(params.aln)

    params.prefix = Channel
        .fromPath(params.aln)
        .map { file -> file.simpleName }

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

    t1_modelfinder_across_rhas_categories(params.prefix, params.aln_ch, params.nthreads)
    keep_modelfinder_results(params.prefix, t1_modelfinder_across_rhas_categories.out[2])
    best_model_per_rhas(params.prefix, keep_modelfinder_results.out)
    best_model_per_rhas.out[2]
        .map { models_out -> store_models(models_out) }
        .set { models_list }
    t1_iqtree_with_best_r2(params.prefix, params.aln, params.nthreads, models_list.map { x -> x[1][0] } )
    get_bic_t1_r2(params.prefix, t1_iqtree_with_best_r2.out[2], "t1_r2")
    hmm_assign_sites_t1_r2(params.prefix, t1_iqtree_with_best_r2.out[5], t1_iqtree_with_best_r2.out[4], "t1_r2")
    evaluate_partitions(params.prefix, hmm_assign_sites_t1_r2.out[1])
    split_aln(params.prefix, params.aln_ch, evaluate_partitions.out[0], params.aln_format)
    t1_iqtree_per_split(split_aln.out[0].flatten(), params.nthreads)
    // TODO: t1_iqtree_class_1 
    // TODO: t1_iqtree_class_2
    // TODO: get_bic_t1_split(.collect())
    concatenate_trees_for_mast(params.prefix, params.aln_ch, t1_iqtree_per_split.out[1].collect(), "class_1_2")
    t2_iqtree_mast(params.prefix, params.aln_ch, params.nthreads, concatenate_trees_for_mast.out[0])
    get_bic_t2_mast(params.prefix, t2_iqtree_mast.out[2], "mast_tr_r2")
    hmm_assign_sites_mast(params.prefix, t2_iqtree_mast.out[4], t2_iqtree_mast.out[6], "mast_tr")
}

