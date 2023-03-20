nextflow.enable.dsl = 2

include {
    test
    update_run_names;
    iqtree_r2;
    hmm_sites_to_ratecats;
    nexus_to_amas;
    evaluate_partitions as evaluate_partitions_1;
    evaluate_partitions as evaluate_partitions_3;
    amas_split as amas_split_1;
    amas_split as amas_split_2;
    iqtree_mfp;
    concat_trees;
    iqtree_hmmster;
    parse_hmmster_partitions
    get_bic;
    store_partitioned_trees;
    prepare_trees;
} from './processes.nf'                            

workflow iteration {
    
    take:
        base_aln
        run_file
        submodel
        nthreads

    main:
        base_aln_ch = Channel.fromPath(base_aln)
        run_file_ch = Channel                     
            .fromPath(run_file)                
            .splitCsv(header:false, sep : "\t")       
            .map{ row -> tuple(row[0], file(row[1])) }

        //test(run_file_ch, base_aln, nthreads)
        split_aln(run_file_ch, base_aln, nthreads)
        //mast(run_file, base_aln, split_aln.out.PartitionedTrees, submodel, nthreads)
        //new_runs = update_run_names(run_file, mast.out.partitioned_aln.collect())
        // TODO: collect BICs and compare
        // TODO: combine NewRuns.tsv

}

workflow split_aln {
    
    /*
     * Split an alignment into two blocks according to R2 FreeRate classes
     */

    take:
        run_file
        base_aln
        nthreads

    main:
        iqtree_r2(run_file, base_aln, nthreads)
        hmm_sites_to_ratecats(run_file, iqtree_r2.out.sitelh, iqtree_r2.out.alninfo)
        nexus_to_amas(run_file, hmm_sites_to_ratecats.out.partitions)
        evaluate_partitions_1(run_file, nexus_to_amas.out.amas_parts)
        amas_split_1(run_file, base_aln, nexus_to_amas.out.amas_parts)
        iqtree_mfp(run_file, amas_split_1.out.aln.flatten(), nthreads)
        store_partitioned_trees(run_file, iqtree_mfp.out.trees.collect())

    emit:
        PartitionedTrees = store_partitioned_trees.out.json

}

workflow mast {
    
    /*
     * Collect all input trees, run HMMSTER, record BIC, make new parts
     */

    take:
        run_name
        base_aln
        partitioned_trees_json
        mast_submodel 
        nthreads

    main:
        prepare_trees(run_name, partitioned_trees_json)
        iqtree_hmmster(run_name, base_aln, prepare_trees.out.input_trees, mast_submodel, nthreads)
        parse_hmmster_partitions(run_name, iqtree_hmmster.out.hmm)
        evaluate_partitions_2(run_name, parse_hmmster_partitions.out.amas_parts)
        amas_split_2(run_name, base_aln, parse_hmmster_partitions.out.amas_parts)

    emit:
        // emit as tuple with run_name
        partitioned_aln = amas_split_2.out.aln.flatten()

}
