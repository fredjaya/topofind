nextflow.enable.dsl = 2

process t1_rk_modelfinder {
    // t1: single-tree
    // rk: use modelfinder to calculate fit of substitution models and k FreeRate categories 
    publishDir "${params.out}/${aln.simpleName}", mode: "copy"

    input:
        path aln
        val nthreads

    output:
        path '*.model.gz'
        path '*.treefile'
        path '*.log'
        path '*.iqtree'
        path '*.ckp.gz'

    script:
    """
    iqtree2 -s ${aln} \
        -pre t1_rk_mf \
        -mrate E,I,R,I+R -m MF\
        -nt ${nthreads}
    """

}

process keep_modelfinder {
    // Keep only ModelFinder results from .log

    publishDir "${params.out}/${aln.simpleName}", mode: "copy"

    input:
        path aln
        path log

    output:
        path models_only

    shell:
    '''
    (sed '0,/^ModelFinder will test up to /d; /^Akaike Information Criterion:/,$d;' | tr -s ' ') < !{log} > models_only
    '''

}

process compare_models {
    // Get overall best fitting model across all Rk
    // Get best fitting model for each Rk

    publishDir "${params.out}/${aln.simpleName}", mode: "copy"

    input:
        path aln
        path models_only

    output:
        path "best_bic_per_rk.tsv"
        path "models_parsed.tsv"
        env MODELS

    script:
    """
    compare_models.R ${models_only}
    MODELS=`sed 1d best_bic_per_rk.tsv`
    """ 
}

def store_models(x) { 
    // Store number of FreeRate categories, model, and BIC scores
    x.tokenize(" ").collate(3) 
}

workflow {

    params.aln = "/home/fredjaya/GitHub/rec-mast/data/data1a.fasta"
    params.out = "/home/fredjaya/Dropbox/treemix_rc/04_testing/mast_sim/"
    params.nthreads = 1
    params.ncpus = 1

    t1_rk_modelfinder(params.aln, params.nthreads)
    keep_modelfinder(params.aln, t1_rk_modelfinder.out[2])
    compare_models(params.aln, keep_modelfinder.out)
    compare_models.out[2]
        .map { models_out -> store_models(models_out) }
        .set { models_list }

    models_list.view()
    // Print the best BICs for each 


    // Proceed if BIC(R2) < BIC(R1)
    //if (compare_models.out[3]) { 
    //    println "R2 BIC is better than R1. Proceed with pipeline :)"
    //}
    //else { 
    //    println "R2 BIC is poorer than R1. Terminiating." 
    //    exit 0
    //}

}

