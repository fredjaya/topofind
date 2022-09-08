nextflow.enable.dsl = 2

process t1_rk_modelfinder {
    // t1: single-tree
    // rk: use modelfinder to calculate fit of substitution model and k FreeRate categories 
    publishDir "${params.out}/${aln.simpleName}"

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

process parse_iqtree_models {
    // Keep only ModelFinder results from .iqtree

    publishDir "${params.out}/${aln.simpleName}"

    input:
        path aln
        path iqtree

    output:
        path models_only

    shell:
    '''
    (sed '0,/^List of models sorted by BIC scores: $/d' | sed '/^AIC, w-AIC   : Akaike information criterion scores and weights.$/,$d' | sed '/^$/d') < !{iqtree} > models_only
    '''

}

workflow {

    params.aln = "/home/fredjaya/GitHub/rec-mast/data/data1a.fasta"
    params.out = "/home/fredjaya/Dropbox/treemix_rc/04_testing/mast_sim/"
    params.nthreads = 1
    params.ncpus = 1

    t1_rk_modelfinder(params.aln, params.nthreads)
    parse_iqtree_models(params.aln, t1_rk_modelfinder.out[3])
}

