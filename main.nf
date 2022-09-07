nextflow.enable.dsl = 2

params.aln = "path_to_aln.fa"

process t1_modelfinder {

    input:
        path aln from params.aln

    output:
        path '*.alninfo' //into run_hmm
        path '*.bionj'
        path '*.ckp.gz'
        path '*.iqtree' //into parse_bic
        path '*.log'
        path '*.mldist'
        path '*.model.gz'
        path '*.sitelh' //into run_hmm
        path '*.siteprob' //into run_hmm
        path '*.treefile'

    script:
    """
    iqtree2 -s ${aln} \
        -pre ${aln.simpleName}/t1_modelfinder \
        -alninfo -wslr -wspr -nt AUTO
    """

}

process t1_r1 {

    input:
        path aln from params.aln

    output:
        path '*.alninfo' //into run_hmm
        path '*.bionj'
        path '*.ckp.gz'
        path '*.iqtree' //into parse_bic
        path '*.log'
        path '*.mldist'
        path '*.model.gz'
        path '*.sitelh' //into run_hmm
        path '*.siteprob' //into run_hmm
        path '*.treefile'

    script:
    """
    iqtree2 -s ${aln} \
        -pre ${aln.simpleName}/t1_r1 \
        -alninfo -wslr -wspr -nt AUTO
    """

}

process t1_r2 {

    input:
        path aln from params.aln

    output:
        path '*.alninfo' //into run_hmm
        path '*.bionj'
        path '*.ckp.gz'
        path '*.iqtree' //into parse_bic
        path '*.log'
        path '*.mldist'
        path '*.model.gz'
        path '*.sitelh' //into run_hmm
        path '*.siteprob' //into run_hmm
        path '*.treefile'

    script:
    """
    iqtree2 -s ${aln} \
        -pre ${aln.simpleName}/t1_r2 \
        -alninfo -wslr -wspr -nt AUTO
    """

}

workflow {

}

