nextflow.enable.dsl = 2

process echo_test {
    input: val run_names
    output: stdout
    script:
    """
    echo "${run_names}"
    """

}
process update_run_names {
    
    input: 
        val run_names

    output: 
        stdout

    script:
    """
    update_run_names.py "${run_names}"
    """

}

process iqtree_r2 {

    publishDir "${params.out}/${run_name}", mode: "copy"

    input:
        val run_name
        path aln
        val nthreads

    output:
        path '*.alninfo', emit: alninfo
        path '*.model.gz'
        path '*.mldist'
        path '*.bionj'
        path '*.sitelh', emit: sitelh
        path '*.treefile', emit: tree
        path '*.siteprob'
        path '*.iqtree'
        path '*.ckp.gz'
        path '*.log'

    script:
    """
    iqtree2 -s ${aln} -pre r2 -mrate E,R2,I+R2 -nt ${nthreads} \
        -wslr -wspr -alninfo  
    """ 
}

process hmm_sites_to_ratecats {
    // Assign sites to FreeRate classes or tree mixtures with the HMM

    publishDir "${params.out}/${run_name}", mode: "copy"

    input:
        val run_name
        path sitelh 
        path alninfo

    output:
        path "site_assignment.png"
        path "r2.partition", emit: partitions

    script:
    """
    hmm_assign_sites.R ${sitelh} ${alninfo} ${run_name}
    """ 
}

process nexus_to_amas {

    publishDir "${params.out}/${run_name}", mode: "copy"

    input:
        val run_name
        path partition 

    output:
        path 'r2.partition_amas', emit: amas_parts

    shell:
    '''
    sed '1,2d; $d; s/\tcharset //; s/;$//' !{partition} > !{partition}_amas
    '''

}

process evaluate_partitions {
    /* 
     * Terminate pipeline if all sites were assigned to a single class
     * i.e. one partition
     * Otherwise, convert .nex partition files output by IQTREE to be 
     * AMAS compliant
     */

    publishDir "${params.out}/${run_name}", mode: "copy"
    debug true
    
    input:
        val run_name
        path partition

    shell:
    '''
    if ((`cat !{partition} | wc -l` == 4)); then
        echo "\nWARNING: All sites in !{run_name} were assigned to a single class."
    fi
    '''
}

process amas_split {
    // Split alignment according to class partitions
   
    publishDir "${params.out}/${run_name}", mode: "copy"

    input:
        val run_name
        path aln
        path partition

    output:
        path "*.fas", emit: aln

    script:
    """
    AMAS.py split -l ${partition} -i ${aln} -f fasta -d dna
    """
}

process iqtree_mfp {

    //errorStrategy { task.exitStatus == 2 ? "ignore" : "terminate" }
    publishDir "${params.out}/${run_name}", mode: "copy"

    input:
        val run_name
        each path(splitted_aln)
        val nthreads

    output:
        path '*.model.gz'
        path '*.treefile', emit: trees
        path '*.log'
        path '*.iqtree'

    script:
    """
    iqtree2 -s ${splitted_aln} -pre ${splitted_aln.simpleName} -nt ${nthreads}
    """
    
}

process concat_trees {
    // Combine trees from partitioned sites for MAST input

    debug true
    publishDir "${params.out}/${run_name}", mode: "copy"
    //errorStrategy "ignore"

    input:
        val run_name
        path aln
        path "*.treefile"

    output:
        path "concat.treefile", optional: true, emit: trees 

    shell:
    '''
    cat *.treefile > concat.treefile
    echo "Concatenating trees for !{run_name}"
    echo "`cat concat.treefile | wc -l` trees will be input to MAST"
    '''
}

process iqtree_hmmster {

    publishDir "${params.out}/${run_name}", mode: "copy"
    errorStrategy { task.exitStatus == 2 ? 'ignore' : 'terminate' } 
    /*
     * exitStatus == 2  
     * ERROR: The number of submodels specified in the mixture does not match 
     * with the tree number 
     * 
     * Due to concat trees having less trees with errorStrategy ignore 
     *
     */

    input:
        val run_name
        path aln
        path trees
        val mast_submodel
        val nthreads

    output:
        path '*.treefile'
        path '*.iqtree', emit: iqtree
        path '*.hmm', emit: hmm
        path '*.ckp.gz'
        path '*.log'

    script:
    """
    iqtree2 -s ${aln} -pre hmmster -te ${trees} -hmmster{gm} \
        -m "TMIX{"${mast_submodel}"}+T" -nt ${nthreads}
    """
    
}

process parse_hmmster_partitions {

    publishDir "${params.out}/${run_name}", mode: "copy"

    input:
        val run_name
        path hmm

    output:
        path 'hmmster.partitions_amas', emit: amas_parts

    script:
    """
    hmmster_partitions.py ${hmm}
    """

}

process get_bic {
    // From *.iqtree for single tree reconstruction and MAST
    // Then append to existing models
    
    publishDir "${params.out}/${run_name}", mode: "copy"

    input:
        val run_name
        path iqtree 

    output:
        env LINE

    shell:
    '''
    MODEL=!{run_name}
    BIC=`grep BIC !{iqtree} | sed -E 's/^.+ //g'`
    LINE=`echo -e "${MODEL}\t${BIC}"`
    '''
}

process compare_bic {

    publishDir "${params.out}/${run_name}", mode: "copy"
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

process store_partitioned_trees {

    publishDir "${params.out}", mode: "copy"
    debug "true"

    input:
        val run_name
        path trees

    output:
        path "PartitionedTrees.json", emit: json
    
    script:
    """
    store_splits.py ${run_name} ${trees}
    """

}

process prepare_trees {

    publishDir "${params.out}/${run_name}", mode: "copy"

    input:
        val run_name
        path PartitionedTrees

    output:
        path "*.tre", emit: input_trees
    
    script:
    """
    retrieve_pt.py ${PartitionedTrees} ${run_name}
    """

}

process bic {

    publishDir "${params.out}", mode: "copy"
    debug "true"

    input:
        val run_name
        path hmm

    output:
        path "BIC.json", emit: json
    
    script:
    """
    store_bic.py ${run_name} ${hmm}
    """

}
