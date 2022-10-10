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
    (sed -E '0,/^ModelFinder will test up to /d; /^Akaike Information Criterion:/,$d; /^WARNING: .+$/d' | tr -s ' ') < !{log} > models_only

    '''

}

process compare_models {
    // Get overall best fitting model across all Rk
    // Get best fitting model for each Rk

    debug "true"
    publishDir "${params.out}/${aln.simpleName}", mode: "copy"

    input:
        path aln
        path models_only

    output:
        path "models_summary.tsv"
        path "modelfinder_parsed.tsv"
        env MODELS

    script:
    """
    compare_models.R ${models_only}
    MODELS=`sed 1d models_summary.tsv`
    cat models_summary.tsv
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
    /* 
     * Terminate pipeline if all sites were assigned to a single class
     * i.e. one partition
     * Otherwise, convert .nex partition files output by IQTREE to be 
     * AMAS compliant
     */

    debug true
    publishDir "${params.out}/${aln.simpleName}", mode: "copy"

    input:
        path aln
        path partition

    output:
        path "*.partition_amas", optional: true

    shell:
    '''
    if ((`cat !{partition} | wc -l` == 4)); then
        echo "\nAll sites in !{aln} were assigned to a single class."
        exit 0
    else
        sed '1,2d; $d; s/\tcharset //; s/;$//' !{partition} > !{partition}_amas
    fi
    '''
}

process split_aln {
    // Split alignment according to class partitions
    
    publishDir "${params.out}/${aln.simpleName}", mode: "copy"

    input:
        path aln
        path partition
        val aln_format

    output:
        path "*.fas"

    script:
    """
    AMAS.py split -l ${partition} -i ${aln} -f ${aln_format} -d dna
    """
}

process iqtree_default {

    publishDir "${params.out}/${aln.simpleName}", mode: "copy"
    errorStrategy { task.exitStatus == 2 ? "ignore" : "terminate" }

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

    debug true
    publishDir "${params.out}/${aln.simpleName}", mode: "copy"

    input:
        path aln
        path "*.treefile"
        val prefix

    output:
        path "*.treefile", optional: true

    shell:
    '''
    #if [ ! -f "*.treefile" ]; then
    #    echo "\nOnly one tree was constructed for !{aln}. Are there enough informative sites in each partition?"
    #    exit 0
    #else
    #    cat *.treefile > !{aln.simpleName}_!{prefix}.treefile
    #fi
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

process get_bic {
    // From *.iqtree for single tree reconstruction and MAST
    // Then append to existing models
    
    debug "true"
    publishDir "${params.out}/${aln.simpleName}", mode: "copy"

    input:
        path aln
        path iqtree 
        val model

    output:
        env LINE

    shell:
    '''
    MODEL=!{model}
    BIC=`grep BIC !{iqtree} | sed -E 's/^.+ //g'`
    LINE=`echo -e "${MODEL}\t${BIC}"`
    echo ${LINE}
    '''
}
