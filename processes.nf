nextflow.enable.dsl = 2

process test {

    debug 'true'
    input:
        tuple val(run_name), file(part_aln)
        path base_aln
        val nthreads
    
    script:
    """
    echo $run_name $part_aln $base_aln $nthreads
    """
}

process iqtree_r2 {

    publishDir "${params.out_dir}/${run_name}", mode: "copy"

    input:
        tuple val(run_name), file(part_aln)
        path base_aln
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
    def aln = part_aln.name != "true_none" ? "${part_aln}" : "${base_aln}"
    """
    iqtree2 -s ${aln} -pre r2 -mrate E+R2,R2,I+R2 -nt ${nthreads} \
        -wslr -wspr -alninfo  
    """ 
}

process hmm_sites_to_ratecats {
    // Assign sites to FreeRate classes or tree mixtures with the HMM

    publishDir "${params.out_dir}/${run_name}", mode: "copy"

    input:
        tuple val(run_name), file(part_aln)
        path sitelh 
        path alninfo

    output:
        path "site_assignment.png"
        path "r2.partition", emit: partitions

    script:
    """
    hmm_assign_sites.R ${sitelh} ${alninfo}
    """ 
}

process nexus_to_amas {

    publishDir "${params.out_dir}/${run_name}", mode: "copy"
    errorStrategy "ignore"
    debug true

    input:
        tuple val(run_name), file(part_aln)
        path partition 

    output:
        path 'r2.partition_amas', emit: amas_parts

    shell:
    '''
    if ((`cat !{partition} | wc -l` == 4)); then
        echo "\nWARNING: All sites in !{run_name} were assigned to a single class."
        exit 1
    else
        sed '1,2d; $d; s/\tcharset //; s/;$//' !{partition} > !{partition}_amas
    fi
    '''

}

process evaluate_partitions {
    /* 
     * Terminate pipeline if all sites were assigned to a single class
     * i.e. one partition
     */

    publishDir "${params.out_dir}/${run_name}", mode: "copy"
    errorStrategy 'ignore'
    debug true
    
    input:
        tuple val(run_name), file(part_aln)
        path partition

    shell:
    '''
    '''
}

process amas_split {
    // Split alignment according to class partitions
   
    publishDir "${params.out_dir}/${run_name}", mode: "copy"

    input:
        tuple val(run_name), file(part_aln)
        path aln
        path partition

    output:
        path "*.fas", emit: aln

    script:
    def aln = part_aln.name != "true_none" ? "${part_aln}" : "${base_aln}"
    """
    AMAS.py split -l ${partition} -i ${aln} -f fasta -d dna
    """
}

process iqtree_mfp {

    //errorStrategy { task.exitStatus == 2 ? "ignore" : "terminate" }
    publishDir "${params.out_dir}/${run_name}", mode: "copy"

    input:
        tuple val(run_name), file(part_aln), file(splitted_aln)
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
    publishDir "${params.out_dir}/${run_name}", mode: "copy"
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

    publishDir "${params.out_dir}/${run_name}", mode: "copy"
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

    publishDir "${params.out_dir}/${run_name}", mode: "copy"

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
    
    publishDir "${params.out_dir}/${run_name}", mode: "copy"

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

    publishDir "${params.out_dir}/${run_name}", mode: "copy"
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

    publishDir "${params.out_dir}", mode: "copy"
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

    publishDir "${params.out_dir}/${run_name}", mode: "copy"

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

    publishDir "${params.out_dir}", mode: "copy"
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

process update_run_names {
   
    publishDir "${params.out_dir}", mode: "copy"
    
    input: 
        val run_name
        val aln_paths

    output: 
        path "NewRuns.tsv"

    script:
    """
    update_run_names.py ${run_name} "${aln_paths}"
    """

}

