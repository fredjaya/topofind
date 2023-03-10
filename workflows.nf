nextflow.enable.dsl = 2

include {                                          
    iqtree_r2;
    hmm_sites_to_ratecats;
    nexus_to_amas;
    evaluate_partitions as evaluate_partitions_1;
    evaluate_partitions as evaluate_partitions_2;
    amas_split as amas_split_1;
    amas_split as amas_split_2;
    iqtree_mfp;
    concat_trees;
    iqtree_hmmster;
    parse_hmmster_partitions
    get_bic;
    store_partitioned_trees;

} from './processes.nf'                            

workflow pre {
    
    take:
        iteration

    main:
        kk
}

workflow split_aln {
    
    /*
     * Split an alignment into two blocks according to R2 FreeRate classes
     */

    take:
        run_name
        aln_ch
        nthreads

    main:
        iqtree_r2(run_name, aln_ch, nthreads)
        hmm_sites_to_ratecats(run_name, iqtree_r2.out.sitelh, iqtree_r2.out.alninfo)
        nexus_to_amas(run_name, hmm_sites_to_ratecats.out.partitions)
        evaluate_partitions_1(run_name, nexus_to_amas.out.amas_parts)
        amas_split_1(run_name, aln_ch, nexus_to_amas.out.amas_parts)
        iqtree_mfp(run_name, amas_split_1.out.aln.flatten(), nthreads)
        store_partitioned_trees(run_name, iqtree_mfp.out.trees.collect())

    emit:
        PartitionedTrees = store_partitioned_trees.out.json

}

workflow mast {
    
    /*
     * Collect all input trees, run HMMSTER, record BIC, make new parts
     */

    take:
        run_name
        aln_ch
        partitioned_trees_json
        mast_submodel 
        nthreads

    main:
        iqtree_hmmster(run_name, aln_ch, trees, mast_submodel, nthreads)
        parse_hmmster_partitions(run_name, iqtree_hmmster.out.hmm)
        evaluate_partitions_2(run_name, parse_hmmster_partitions.out.amas_parts)
        amas_split_2(run_name, aln_ch, parse_hmmster_partitions.out.amas_parts)
        //amas_split.out[0]
        //    .flatten()
        //    .branch {
        //        class_1: it =~ /class_1/
        //        class_2: it =~ /class_2/
        //        class_3: it =~ /class_3/
        //    } .set { splitted_aln }

    //emit:
    //    class_1 = splitted_aln.class_1
    //    class_2 = splitted_aln.class_2
    //    class_3 = splitted_aln.class_3
    //    bic = get_bic_t2_mast.out[0]

}

workflow evaluate_run {

    take:    

    main:
        get_bic(run_name, iqtree_hmmster.out.iqtree)

}
