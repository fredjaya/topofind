nextflow.enable.dsl = 2

params.aln_path = "/home/fredjaya/Dropbox/treemix_rc/04_testing/rob_sims/test1.fa"
params.aln_format = "fasta"
params.out = "/home/fredjaya/Dropbox/treemix_rc/04_testing/rob_sims/"
params.nthreads = 1
params.ncpus = 1

params.aln = Channel.fromPath(params.aln_path)

def store_models(x) {
    // Store number of FreeRate categories, model, and BIC scores
    x.tokenize(" ").collate(3)
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
    t1_rk_modelfinder; 
    keep_modelfinder; compare_models; 
    t1_r2;
    hmm_assign_sites; 
    parse_partition; 
    split_aln; 
    iqtree_default;
    concatenate_trees; 
    iqtree_mast_tr_r2; 
    hmm_assign_sites as hmm_assign_sites_mast_tr
} from './processes.nf' 

workflow {
    
    /*
     * This workflow will only use R2 rate categories to to delimit sites 
     */ 
    
    log.info"""
    base        = ${baseDir}
    out         = ${params.out}
    nthreads    = ${params.nthreads}
    ncpus       = ${params.ncpus}
    aln_path    = ${params.aln_path}
    aln_format  = ${params.aln_format}
    """

    t1_rk_modelfinder(params.aln, params.nthreads)
    keep_modelfinder(params.aln, t1_rk_modelfinder.out[2])
    compare_models(params.aln, keep_modelfinder.out)
    compare_models.out[2]
        .map { models_out -> store_models(models_out) }
        .set { models_list }

    println("\nBest model per number of FreeRates categories:")
    models_list.view()

    /*
    // Compares BICs
    models_list
        .map { models_list -> compare_bic(models_list, 0, 1) }
        .view()
    */

    // For now, just proceed with R2

    t1_r2(params.aln, params.nthreads, models_list.map { x -> x[1][1] } )
    hmm_assign_sites(params.aln, t1_r2.out[5], t1_r2.out[4], "t1_r2")
    parse_partition(params.aln, hmm_assign_sites.out[1])
    split_aln(params.aln, parse_partition.out[0], params.aln_format)
    iqtree_default(split_aln.out[0].flatten(), params.nthreads)
    concatenate_trees(params.aln, iqtree_default.out[1].collect(), "class_1_2")
    iqtree_mast_tr_r2(params.aln, params.nthreads, concatenate_trees.out[0])
    hmm_assign_sites_mast_tr(params.aln, iqtree_mast_tr_r2.out[4], iqtree_mast_tr_r2.out[6], "mast_tr")
}

