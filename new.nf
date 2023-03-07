nextflow.enable.dsl = 2

include {

    split_aln;

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

    split_aln(params.aln_name, "01_split_aln", params.aln_ch, params.nthreads)

}
