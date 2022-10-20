nextflow.enable.dsl = 2

include {                                          
    t1_modelfinder_across_rhas_categories;         
    keep_modelfinder_results;
    best_model_per_rhas; 
    store_models;
    t1_iqtree_with_best_r2;                        
    get_bic as get_bic_t1_r2;                      
    hmm_assign_sites as hmm_assign_sites_t1_r2;    
    evaluate_partitions;                           
    amas_split;                                     
    t1_iqtree_per_split;                           
    get_bic as get_bic_t1_split;                   
    concatenate_trees_for_mast;                    
    t2_iqtree_mast;                                
    get_bic as get_bic_t2_mast;                    
    hmm_assign_sites as hmm_assign_sites_mast;     
} from './processes.nf'                            

workflow split_aln {
    
    /*
     * Split an alignment into two blocks according to R2 FreeRate classes
     *
     */
    take:
        aln_name
        run_mode
        aln_ch
        aln_format
        nthreads

    main:
        t1_modelfinder_across_rhas_categories(aln_name, run_mode, aln_ch, nthreads)
        keep_modelfinder_results(aln_name, run_mode, t1_modelfinder_across_rhas_categories.out[2])
        best_model_per_rhas(aln_name, run_mode, keep_modelfinder_results.out)
        best_model_per_rhas.out[2]
            .map { models_out -> store_models(models_out) }
            .set { models_list }
        t1_iqtree_with_best_r2(aln_name, run_mode, aln_ch, nthreads, models_list.map { x -> x[1][0] } )
        get_bic_t1_r2(aln_name, run_mode, t1_iqtree_with_best_r2.out[2])
        hmm_assign_sites_t1_r2(aln_name, run_mode, t1_iqtree_with_best_r2.out[5], t1_iqtree_with_best_r2.out[4])
        evaluate_partitions(aln_name, run_mode, hmm_assign_sites_t1_r2.out[1])
        amas_split(aln_name, run_mode, aln_ch, evaluate_partitions.out[0], aln_format)
        t1_iqtree_per_split(aln_name, run_mode, amas_split.out[0].flatten(), nthreads)

    emit:
        t1 = t1_iqtree_with_best_r2.out[0]
        t2 = t1_iqtree_per_split.out[1].collect()
        bic = get_bic_t1_r2.out[0]
}

workflow mast {
    
    /*
     * Collect all input trees, run MAST +TR, record the BIC and 
     * assign sites with the HMM
     */

    take:
        aln_name
        run_mode
        aln_ch
        aln_format
        trees
        nthreads

    main:
        concatenate_trees_for_mast(aln_name, run_mode, aln_ch, trees)
        t2_iqtree_mast(aln_name, run_mode, aln_ch, nthreads, concatenate_trees_for_mast.out[0])
        get_bic_t2_mast(aln_name, run_mode, t2_iqtree_mast.out[2])
        hmm_assign_sites_mast(aln_name, run_mode, t2_iqtree_mast.out[4], t2_iqtree_mast.out[6])
        evaluate_partitions(aln_name, run_mode, hmm_assign_sites_mast.out[1])
        amas_split(aln_name, run_mode, aln_ch, evaluate_partitions.out[0], aln_format)
        amas_split.out[0]
            .flatten()
            .branch {
                class_1: it =~ /class_1/
                class_2: it =~ /class_2/
            } .set { splitted_aln }

    emit:
        class_1 = splitted_aln.class_1
        class_2 = splitted_aln.class_2
        bic = get_bic_t2_mast.out[0]

}

workflow sort_trees {

    /*
     * Collect t1 and split t2 trees. Branch channels to contain three trees.
     * i.e only one split (three blocks) is analysed at once.
     */

    take:
        mast_tree_class_1
        mast_tree_class_2
        new_trees

    main:
        new_trees
            .branch {
                class_1_split: it =~ /class_1-out_class/
                class_2_split: it =~ /class_2-out_class/
            } .set { sorted_trees }

        mast_tree_class_1.mix(sorted_trees.class_2_split).set { trees_a1b12 }
        mast_tree_class_2.mix(sorted_trees.class_1_split).set { trees_a12b1 }

    emit:
        trees_a1b12
        trees_a12b1

}
