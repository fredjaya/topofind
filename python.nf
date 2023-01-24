nextflow.enable.dsl = 2

include {

    split_aln;
    mast;

} from "./workflows.nf"

params.aln_ch = Channel
    .fromPath(params.aln)

params.aln_name = Channel
    .fromPath(params.aln)
    .map { file -> file.simpleName }

workflow {

    if( params.mode == "first_iter" ) {
        split_aln(params.aln_name, params.run_name, params.aln_ch, params.aln_format, params.nthreads)
        mast(params.aln_name, params.run_name, params.aln_ch, params.aln_format, split_aln.out.t2, "GTR+FO+G,GTR+FO+G", params.nthreads)
    }

    if( params.mode == "split_aln" ) {        
        params.trees_ch = Channel
            .fromPath(params.trees)
        split_aln(params.aln_name, params.run_name, params.aln_ch, params.aln_format, params.nthreads)
    }

    if( params.mode == "mast" ) {
        params.trees_ch = Channel
            .fromPath(params.trees)
        mast(params.aln_name, params.run_name, params.aln_ch, params.aln_format, params.trees_ch, params.submodel, params.nthreads)
    }
}
