process BUILD_FASTA_INDEX {
    tag 'reference-faidx'
    publishDir "${params.outdir}/reference", mode: 'copy', pattern: 'reference.fa*'

    input:
    path(reference_source)

    output:
    tuple path('reference.fa'), path('reference.fa.fai'), emit: reference

    script:
    """
    set -euo pipefail
    cat ${reference_source} > reference.fa
    samtools faidx reference.fa
    """

    stub:
    """
    cat ${reference_source} > reference.fa
    touch reference.fa.fai
    """
}

process BUILD_SEQUENCE_DICTIONARY {
    tag 'reference-dict'
    publishDir "${params.outdir}/reference", mode: 'copy', pattern: 'reference.dict'

    input:
    path(reference_source)

    output:
    path('reference.dict'), emit: dict

    script:
    """
    set -euo pipefail
    cat ${reference_source} > reference.fa
    gatk CreateSequenceDictionary -R reference.fa -O reference.dict
    """

    stub:
    """
    touch reference.dict
    """
}

process BUILD_BWAMEM2_INDEX {
    tag 'reference-bwamem2'
    publishDir "${params.outdir}/reference/bwamem2", mode: 'copy'

    input:
    path(reference_source)

    output:
    tuple path('reference.fa'), path('reference.fa.*'), emit: reference

    script:
    """
    set -euo pipefail
    cat ${reference_source} > reference.fa
    bwa-mem2 index reference.fa
    """

    stub:
    """
    cat ${reference_source} > reference.fa
    touch reference.fa.0123 reference.fa.amb reference.fa.ann reference.fa.bwt.2bit.64 reference.fa.pac
    """
}