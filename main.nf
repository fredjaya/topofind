nextflow.enable.dsl = 2

process t1_r1 {
    // Single-tree
    // Estimate best substitution model without RHAS
    publishDir "${params.out}/${aln.simpleName}"

    input:
        path aln
        val nthreads

    output:
        path '*.alninfo'
        path '*.bionj'
        path '*.ckp.gz'
        path '*.iqtree'
        path '*.log'
        path '*.mldist'
        path '*.model.gz'
        path '*.sitelh'
        path '*.treefile'

    script:
    """
    iqtree2 -s ${aln} \
        -pre t1_r1 \
        -mrate E,F,I \
        -alninfo -wslr -nt ${nthreads}
    """

}

process t1_r2 {

    input:
        path aln from params.aln

    output:
        path '*.alninfo'
        path '*.bionj'
        path '*.ckp.gz'
        path '*.iqtree'
        path '*.log'
        path '*.mldist'
        path '*.model.gz'
        path '*.sitelh'
        path '*.siteprob'
        path '*.treefile'

    script:
    """
    iqtree2 -s ${aln} \
        -pre t1_r2 \
        -alninfo -wslr -wspr -nt ${nthreads}
    """

}

workflow {

    params.aln = "/home/fredjaya/Dropbox/treemix_rc/04_testing/mast_sim/data1a.fasta"
    params.out = "/home/fredjaya/Dropbox/treemix_rc/04_testing/mast_sim/"
    params.nthreads = 1
    params.ncpus = 1

    t1_r1(params.aln, params.nthreads)
    println "${params.out}"

}

