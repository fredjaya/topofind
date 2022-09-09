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
        path models_only

    output:
        path models_parsed
        path best_bic_rk

    script:
    """
    Rscript parse_models.R ${models_only}
    """
}

workflow {

    params.aln = "/home/fredjaya/GitHub/rec-mast/data/data1a.fasta"
    params.out = "/home/fredjaya/Dropbox/treemix_rc/04_testing/mast_sim/"
    params.nthreads = 1
    params.ncpus = 1

    t1_rk_modelfinder(params.aln, params.nthreads)
    keep_modelfinder(params.aln, t1_rk_modelfinder.out[2])
    //parse_models()

    // Proceed if BIC(R2) < BIC(R1)
}

