nextflow.enable.dsl = 2

include {                                          
    iqtree_r2;
    hmm_sites_to_ratecats;
    evaluate_partitions;
    amas_split as amas_split_1;
    amas_split as amas_split_2;
    iqtree_mfp;
    concat_trees;
    iqtree_hmmster;
    parse_hmmster_partitions
    get_bic;
} from './processes.nf'                            


workflow split_aln {
    
    /*
     * Split an alignment into two blocks according to R2 FreeRate classes
     */

    take:
        aln_name
        run_name
        aln_ch
        nthreads

    main:
        iqtree_r2(aln_name, run_name, aln_ch, nthreads)
        hmm_sites_to_ratecats(aln_name, run_name, iqtree_r2.out.sitelh, iqtree_r2.out.alninfo)
        evaluate_partitions(aln_name, run_name, hmm_sites_to_ratecats.out.partitions)
        amas_split_1(aln_name, run_name, aln_ch, evaluate_partitions.out.amas_parts)
        iqtree_mfp(aln_name, run_name, amas_split_1.out[0].flatten(), nthreads)

    emit:
        t1 = iqtree_r2.out.tree
        t2 = iqtree_mfp.out.trees.collect()
}

workflow mast {
    
    /*
     * Collect all input trees, run HMMSTER, record BIC, make new parts
     */

    take:
        aln_name
        run_name
        aln_ch
        trees
        mast_submodel 
        nthreads

    main:
        concat_trees(aln_name, run_name, aln_ch, trees)
        iqtree_hmmster(aln_name, run_name, aln_ch, concat_trees.out.trees, mast_submodel, nthreads)
        parse_hmmster_partitions(aln_name, iqtree_hmmster.out.hmm)
        //evaluate_partitions(aln_name, run_name, iqtree.out[1])
        //amas_split_2(aln_name, run_name, aln_ch, evaluate_partitions.out[0], aln_format)
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
        get_bic(aln_name, run_name, iqtree_hmmster.out.iqtree)

}
