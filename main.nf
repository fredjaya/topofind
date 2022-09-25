nextflow.enable.dsl = 2

process t1_rk_modelfinder {
    /*
     * t1: single-tree
     * rk: use modelfinder to calculate fit of substitution models and 
     *     k FreeRate categories
     */

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

process t1_r2 {

    /*
     * t1: single-tree
     * r2: k=2 rate categories
     */

    publishDir "${params.out}/${aln.simpleName}", mode: "copy"

    input:
        path aln
        val nthreads
        val model

    output:
        path '*.treefile'
        path '*.log'
        path '*.iqtree'
        path '*.ckp.gz'
        path '*.alninfo'
        path '*.sitelh'
        path '*.siteprob'

    script:
    """
    iqtree2 -s ${aln} \
        -pre t1_r2 \
        -m ${model} \
        -wslr -wspr -alninfo \
        -nt ${nthreads}
    """

}

process hmm_assign_sites {
    // Assign sites to FreeRate classes with the HMM

    publishDir "${params.out}/${aln.simpleName}",
        saveAs: { filename -> "t1_r2$filename" },
        mode: "copy"

    input:
        path aln
        path sitelh 
        path alninfo

    output:
        path "_site_assignment.png"
        path ".partition"

    script:
    """
    hmm_assign_sites.R ${sitelh} ${alninfo}
    """ 
}

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

workflow {
    
    /*
     * This workflow will only use R2 rate categories to to delimit sites 
     */ 

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
    hmm_assign_sites(params.aln, t1_r2.out[5], t1_r2.out[4])
}

