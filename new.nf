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

    log.info"""
    base        = ${baseDir}
    out         = ${params.out}
    aln         = ${params.aln}
    nthreads    = ${params.nthreads}
    """

    split_aln("01_split_A_B", params.aln_ch, params.nthreads)
    mast("02_mast_A_B", params.aln_ch, split_aln.out.t2, "GTR+FO+G,GTR+FO+G", params.nthreads)
}
