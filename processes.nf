nextflow.enable.dsl = 2

def store_models(x) {
    // Store number of FreeRate categories, model, and BIC scores
    x.tokenize(" ").collate(2)                                   
}

//def count_trees = 
process t1_modelfinder_across_rhas_categories {
    /*
     * t1: single-tree
     * rk: use modelfinder to calculate fit of substitution models and 
     *     k FreeRate categories
     */
    
    publishDir "${params.out}/${aln_name}/${run_mode}", mode: "copy"

    input:
        val aln_name
        val run_mode
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

process keep_modelfinder_results {
    // Keep only ModelFinder results from .log

    publishDir "${params.out}/${aln_name}/${run_mode}", mode: "copy"

    input:
        val aln_name
        val run_mode
        path log

    output:
        path models_only

    shell:
    '''
    (sed -E '0,/^ModelFinder will test up to /d; /^Akaike Information Criterion:/,$d; /^WARNING: .+$/d' | tr -s ' ') < !{log} > models_only

    '''

}

process best_model_per_rhas {
    // Get overall best fitting model across all Rk
    // Get best fitting model for each Rk

    debug "true"

    publishDir "${params.out}/${aln_name}/${run_mode}", mode: "copy"

    input:
        val aln_name
        val run_mode
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

process t1_iqtree_with_best_r2 {

    /*
     * t1: single-tree
     * r2: k=2 rate categories
     */

    publishDir "${params.out}/${aln_name}/${run_mode}", mode: "copy"

    input:
        val aln_name
        val run_mode
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

    publishDir "${params.out}/${aln_name}/${run_mode}", mode: "copy"

    input:
        val aln_name 
        val run_mode
        path sitelh 
        path alninfo

    output:
        path "*_site_assignment.png"
        path "*.partition"

    script:
    """
    hmm_assign_sites.R ${sitelh} ${alninfo} ${run_mode}
    """ 
}

process evaluate_partitions {
    /* 
     * Terminate pipeline if all sites were assigned to a single class
     * i.e. one partition
     * Otherwise, convert .nex partition files output by IQTREE to be 
     * AMAS compliant
     */

    publishDir "${params.out}/${aln_name}/${run_mode}", mode: "copy"
    debug true
    
    input:
        val aln_name
        val run_mode
        path partition

    output:
        path "*.partition_amas", optional: true

    shell:
    '''
    if ((`cat !{partition} | wc -l` == 4)); then
        echo "\nWARNING: All sites in !{run_mode} were assigned to a single class."
    fi
        sed '1,2d; $d; s/\tcharset //; s/;$//' !{partition} > !{partition}_amas
    '''
}

process amas_split {
    // Split alignment according to class partitions
   
    publishDir "${params.out}/${aln_name}/${run_mode}", mode: "copy"

    input:
        val aln_name
        val run_mode
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

process t1_iqtree_per_split {

    //errorStrategy { task.exitStatus == 2 ? "ignore" : "terminate" }
    publishDir "${params.out}/${aln_name}/${run_mode}", mode: "copy"

    input:
        val aln_name
        val run_mode
        each path(splitted_aln)
        val nthreads

    output:
        path '*.model.gz'
        path '*.treefile'
        path '*.log'
        path '*.iqtree'

    script:
    """
    iqtree2 -s ${splitted_aln} -pre ${splitted_aln.simpleName} -nt ${nthreads}
    """
    
}

process concatenate_trees_for_mast {
    // Combine trees from partitioned sites for MAST input

    debug true
    publishDir "${params.out}/${aln_name}/${run_mode}", mode: "copy"
    errorStrategy "ignore"

    input:
        val aln_name
        val run_mode
        path aln
        path "*.treefile"

    output:
        path "*.treefile", optional: true

    shell:
    '''
    cat *.treefile > concat.treefile
    '''
}

process iqtree_mast_t_r2 {

    publishDir "${params.out}/${aln_name}/${run_mode}", mode: "copy"

    input:
        path aln
        val run_mode
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

process t2_iqtree_mast {

    publishDir "${params.out}/${aln_name}/${run_mode}", mode: "copy"

    input:
        val aln_name
        val run_mode
        path aln
        path trees
        val mast_submodel
        val nthreads

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
    iqtree2 -s ${aln} -pre t2_mast_tr -nt ${nthreads} -te ${trees} \
        -m "TMIX{"${mast_submodel}"}+TR" -nt ${nthreads} -wslr -wspr -alninfo
    """
    
}

process get_bic {
    // From *.iqtree for single tree reconstruction and MAST
    // Then append to existing models
    
    publishDir "${params.out}/${aln_name}/${run_mode}", mode: "copy"

    input:
        val aln_name 
        val run_mode
        path iqtree 

    output:
        env LINE

    shell:
    '''
    MODEL=!{run_mode}
    BIC=`grep BIC !{iqtree} | sed -E 's/^.+ //g'`
    LINE=`echo -e "${MODEL}\t${BIC}"`
    '''
}

process compare_bic {

    publishDir "${params.out}/${aln_name}/${run_mode}", mode: "copy"
    debug "true"

    input:
        tuple val(step_0), val(bic_0)
        tuple val(step_1), val(bic_1)

   
    output:
        env OUT_STEP 

    shell:
    '''
    OUT_STEP=`compare_bic.R !{step_0} !{bic_0} !{step_1} !{bic_1}`
    '''

}
