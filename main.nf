nextflow.enable.dsl = 2

include {

    split_aln as split_aln_1;

    mast as mast_2A_2B;

    split_aln as split_aln_2A;
    split_aln as split_aln_2B;

    mast as mast_2A_3BA_3BB;
    mast as mast_3AA_3AB_2B;

    split_aln as split_aln_3AA; 
    split_aln as split_aln_3AB; 
    split_aln as split_aln_3BA; 
    split_aln as split_aln_3BB; 

    mast as mast_2A_3BA_4BBA_4BBB; 
    mast as mast_2A_4BAA_4BAB_3BB; 
    mast as mast_4AAA_4AAB_3AB_2B; 
    mast as mast_3AA_4ABA_4ABB_2B;
    mast as mast_3AA_3AB_3BA_3BB

} from "./workflows.nf"

include { compare_bic } from "./processes.nf"

params.aln_ch = Channel
    .fromPath(params.aln)

params.aln_name = Channel
    .fromPath(params.aln)
    .map { file -> file.simpleName }

workflow {

    log.info"""
    // Paths and directories
    base        = ${baseDir}
    out         = ${params.out}

    // Resource usage
    nthreads    = ${params.nthreads}

    // Input alignment
    aln         = ${params.aln}
    aln_format  = ${params.aln_format}

    // Debugging
    previous_model = ${params.previous_model}
    """

    if( params.mode = "split_aln" )
	    log.info"""split_aln"""
    else
    	log.info"""

	***********************************
	* RUNNING WHOLE PIPELINE MANUALLY *
	***********************************
	"""

    	/*
    	 * 1-to-2 block split
    	 */
    	split_aln_1(params.aln_name, "01_split_aln", params.aln_ch, params.aln_format, params.nthreads)
    	split_aln_1.out.bic.map{ it -> it.tokenize(" ") }.set{ bic_1 }

    	/*
    	 * 2-tree MAST
    	 */
    	mast_2A_2B(params.aln_name, "02_mast_2A_2B", params.aln_ch, params.aln_format, split_aln_1.out.t2, "GTR+FO+G,GTR+FO+G", params.nthreads)
    	mast_2A_2B.out.bic.map{ it -> it.tokenize(" ") }.set{ bic_2 }

    	/*
    	 * 2_to_3 block splits
    	 */
    	split_aln_2A(params.aln_name, "03_split_aln_2A", mast_2A_2B.out.class_1, "fasta", params.nthreads)
    	split_aln_2B(params.aln_name, "03_split_aln_2B", mast_2A_2B.out.class_2, "fasta", params.nthreads)

    	/*
    	 * 3_tree MAST
    	 */
    	split_aln_2A.out.t1.concat(split_aln_2B.out.t2).collect().set { trees_2A_3BA_3BB }
    	split_aln_2A.out.t2.concat(split_aln_2B.out.t1).collect().set { trees_3AA_3AB_2B }

    	mast_2A_3BA_3BB(params.aln_name, "04_mast_2A_3BA_3BB", params.aln_ch, params.aln_format, trees_2A_3BA_3BB, "GTR+FO+G,GTR+FO+G,GTR+FO+G", params.nthreads)
    	mast_3AA_3AB_2B(params.aln_name, "04_mast_3AA_3AB_2B", params.aln_ch, params.aln_format, trees_3AA_3AB_2B, "GTR+FO+G,GTR+FO+G,GTR+FO+G", params.nthreads)

    	mast_2A_3BA_3BB.out.bic.map { it -> it.tokenize(" ") }.set{ bic_3a }
    	mast_3AA_3AB_2B.out.bic.map { it -> it.tokenize(" ") }.set{ bic_3b }

    	/*
    	 * 3_to_4 block splits
    	 */
    	//2A split == 3AA_3AB
    	split_aln_3AA(params.aln_name, "05_split_aln_3AA", mast_3AA_3AB_2B.out.class_1, "fasta", params.nthreads)
    	split_aln_3AB(params.aln_name, "05_split_aln_3AB", mast_3AA_3AB_2B.out.class_2, "fasta", params.nthreads)
    	//2B split == 3BA_3BB
    	split_aln_3BA(params.aln_name, "05_split_aln_3BA", mast_2A_3BA_3BB.out.class_1, "fasta", params.nthreads)
    	split_aln_3BB(params.aln_name, "05_split_aln_3BB", mast_2A_3BA_3BB.out.class_2, "fasta", params.nthreads)

   /	*
    	* 4_tree MAST
    	*/
    	trees_2A_3BA_4BBA_4BBB =
    	    split_aln_2A.out.t1
    	    .concat(split_aln_3BA.out.t1)
    	    .concat(split_aln_3BB.out.t2)
    	    .collect()
    	trees_2A_4BAA_4BAB_3BB =
    	    split_aln_2A.out.t1
    	    .concat(split_aln_3BA.out.t2)
    	    .concat(split_aln_3BB.out.t1)
    	    .collect()
    	trees_4AAA_4AAB_3AB_2B =
    	    split_aln_3AA.out.t2
    	    .concat(split_aln_3AB.out.t1)
    	    .concat(split_aln_2B.out.t1)
    	    .collect()
    	trees_2A_4BAA_4BAB_3BB =
    	    split_aln_2A.out.t1
    	    .concat(split_aln_3BA.out.t2)
    	    .concat(split_aln_3BB.out.t1)
    	    .collect()
    	trees_3AA_4ABA_4ABB_2B =
    	    split_aln_3AA.out.t1
    	    .concat(split_aln_3AB.out.t2)
    	    .concat(split_aln_2B.out.t1)
    	    .collect()
    	trees_3AA_3AB_3BA_3BB =
    	    split_aln_2A.out.t2
    	    .concat(split_aln_3BA.out.t2)
    	    .collect()

    	mast_2A_3BA_4BBA_4BBB(params.aln_name, "06_mast_2A_3BA_4BBA_4BBB", params.aln_ch, params.aln_format, trees_2A_3BA_4BBA_4BBB, "GTR+FO+G,GTR+FO+G,GTR+FO+G,GTR+FO+G", params.nthreads)
    	mast_2A_4BAA_4BAB_3BB(params.aln_name, "06_mast_2A_4BAA_4BAB_3BB", params.aln_ch, params.aln_format, trees_2A_4BAA_4BAB_3BB, "GTR+FO+G,GTR+FO+G,GTR+FO+G,GTR+FO+G", params.nthreads)
    	mast_4AAA_4AAB_3AB_2B(params.aln_name, "06_mast_4AAA_4AAB_3AB_2B", params.aln_ch, params.aln_format, trees_4AAA_4AAB_3AB_2B, "GTR+FO+G,GTR+FO+G,GTR+FO+G,GTR+FO+G", params.nthreads)
    	mast_3AA_4ABA_4ABB_2B(params.aln_name, "06_mast_3AA_4ABA_4ABB_2B", params.aln_ch, params.aln_format, trees_3AA_4ABA_4ABB_2B, "GTR+FO+G,GTR+FO+G,GTR+FO+G,GTR+FO+G", params.nthreads)
    	mast_3AA_3AB_3BA_3BB(params.aln_name, "06_mast_3AA_3AB_3BA_3BB", params.aln_ch, params.aln_format, trees_3AA_3AB_3BA_3BB, "GTR+FO+G,GTR+FO+G,GTR+FO+G,GTR+FO+G", params.nthreads)

    	/*
    	 * Print BIC scores
    	 */
    	bic_1
    	    .concat(bic_2)
    	    .concat(bic_3a)
    	    .concat(bic_3b)
    	    .collect()
    	    .map { it -> it.collate(2) } 
    	    .set { bic_all }

    	bic_all.view()

}
