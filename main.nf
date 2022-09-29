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
    // Assign sites to FreeRate classes or tree mixtures with the HMM

    publishDir "${params.out}/${aln.simpleName}", mode: "copy"

    input:
        path aln
        path sitelh 
        path alninfo
        val prefix

    output:
        path "*_site_assignment.png"
        path "*.partition"

    script:
    """
    hmm_assign_sites.R ${sitelh} ${alninfo} ${prefix}
    """ 
}

process parse_partition {
    // Convert .nex partition files output by IQTREE to be compatible with AMAS

    publishDir "${params.out}/${aln.simpleName}", mode: "copy"

    input:
        path aln
        path partition

    output:
        path "*.partition_amas"

    shell:
    '''
    echo !{partition}
    sed '1,2d; $d; s/\tcharset //; s/;$//' !{partition} > !{partition}_amas
    '''
}

process split_aln {
    // Split alignment according to class partitions
    
    publishDir "${params.out}/${aln.simpleName}", mode: "copy"

    input:
        path aln
        path partition

    output:
        path "*.fas"

    script:
    """
    AMAS.py split -l ${partition} -i ${aln} -f fasta -d dna
    """
}

process iqtree_default {

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
    iqtree2 -s ${aln} -pre ${aln.simpleName} -nt ${nthreads}
    """
    
}

process concatenate_trees {
    // Combine trees from partitioned sites for MAST input

    publishDir "${params.out}/${aln.simpleName}", mode: "copy"

    input:
        path aln
        path "*.treefile"
        val prefix

    output:
        path "*.treefile"

    shell:
    '''
    cat *.treefile > !{aln.simpleName}_!{prefix}.treefile
    '''
}

process iqtree_mast_t_r2 {

    publishDir "${params.out}/${aln.simpleName}", mode: "copy"

    input:
        path aln
        val nthreads
        path trees

    output:
        path '*.treefile'
        path '*.log'
        path '*.iqtree'
        path '*.ckp.gz'
        path '*.sitelh'
        path '*.siteprob'
        path '*.alninfo'

    script:
    """
    iqtree2 -s ${aln} -pre ${trees.simpleName}_mast_t -nt ${nthreads} -te ${trees} \
        -m "TMIX{GTR+FO+G,GTR+FO+G}+T" -nt ${nthreads} -wslr -wspr -alninfo
    """
    
}

process iqtree_mast_tr_r2 {

    publishDir "${params.out}/${aln.simpleName}", mode: "copy"

    input:
        path aln
        val nthreads
        path trees

    output:
        path '*.treefile'
        path '*.log'
        path '*.iqtree'
        path '*.ckp.gz'
        path '*.sitelh'
        path '*.siteprob'
        path '*.alninfo'

    script:
    """
    iqtree2 -s ${aln} -pre ${trees.simpleName}_mast_tr -nt ${nthreads} -te ${trees} \
        -m "TMIX{GTR+FO+G,GTR+FO+G}+TR" -nt ${nthreads} -wslr -wspr -alninfo
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
    hmm_assign_sites(params.aln, t1_r2.out[5], t1_r2.out[4], "t1_r2")
    parse_partition(params.aln, hmm_assign_sites.out[1])
    split_aln(params.aln, parse_partition.out[0])
    iqtree_default(split_aln.out[0].flatten(), params.nthreads)
    concatenate_trees(params.aln, iqtree_default.out[1].collect(), "class_1_2")
    iqtree_mast_t_r2(params.aln, params.nthreads, concatenate_trees.out[0])
    iqtree_mast_tr_r2(params.aln, params.nthreads, concatenate_trees.out[0])
    //hmm_assign_sites(params.aln, iqtree_mast_t_r2.out[4], iqtree_mast_t_r2.out[6], "mast_t")
    //hmm_assign_sites(params.aln, iqtree_mast_tr_r2.out[4], iqtree_mast_tr_r2.out[6], "mast_tr")
}

