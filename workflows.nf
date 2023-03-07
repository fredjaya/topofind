nextflow.enable.dsl = 2

include {                                          
    iqtree_r2;
    hmm_sites_to_ratecats;
    evaluate_partitions;
    amas_split;
    iqtree_mfp;
} from './processes.nf'                            


workflow split_aln {
    
    /*
     * Split an alignment into two blocks according to R2 FreeRate classes
     *
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
        amas_split(aln_name, run_name, aln_ch, evaluate_partitions.out.amas_parts)
        iqtree_mfp(aln_name, run_name, amas_split.out[0].flatten(), nthreads)

    emit:
        t1 = iqtree_r2.out.tree
        t2 = iqtree_mfp.out.trees.collect()
}

workflow mast {
    
    /*
     * Collect all input trees, run MAST +TR, record the BIC and 
     * assign sites with the HMM
     */

    take:
        aln_name
        run_name
        aln_ch
        aln_format
        trees
        mast_submodel 
        nthreads

    main:
        t2_iqtree_mast(aln_name, run_name, aln_ch, trees, mast_submodel, nthreads)
        get_bic_t2_mast(aln_name, run_name, t2_iqtree_mast.out[2])
        hmm_assign_sites_mast(aln_name, run_name, t2_iqtree_mast.out[4], t2_iqtree_mast.out[6])
        evaluate_partitions(aln_name, run_name, hmm_assign_sites_mast.out[1])
        amas_split(aln_name, run_name, aln_ch, evaluate_partitions.out[0], aln_format)
        amas_split.out[0]
            .flatten()
            .branch {
                class_1: it =~ /class_1/
                class_2: it =~ /class_2/
                class_3: it =~ /class_3/
            } .set { splitted_aln }

    emit:
        class_1 = splitted_aln.class_1
        class_2 = splitted_aln.class_2
        class_3 = splitted_aln.class_3
        bic = get_bic_t2_mast.out[0]

}
