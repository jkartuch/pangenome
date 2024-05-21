//
// Check input FASTA and prepare indices
//

include { GRABIX_CHECK                } from '../../modules/nf-core/grabix/check/main.nf'
include { SAMTOOLS_FAIDX              } from '../../modules/nf-core/samtools/faidx/main.nf'
include { TABIX_BGZIP                 } from '../../modules/nf-core/tabix/bgzip/main.nf'

workflow INPUT_CHECK {
    take:
    fasta // file: /path/to/sequences.fasta

    main:

    ch_versions = Channel.empty() // we collect all versions here
    ch_fasta = Channel.empty() // final output channel [ val(meta) , [ fasta ] ]

    fai_path = file("${params.input}.fai")
    gzi_path = file("${params.input}.gzi")

    fai = Channel.empty() // we store the .fai index here [ fai ]
    gzi = Channel.empty() // we store the .gzi index here [ gzi ]

    meta_fasta = Channel.empty() // intermediate channel where we build our [ val(meta) , [ fasta ] ]
    fasta_file_name = fasta.getName()


    meta_fasta = tuple([ id:fasta_file_name ], fasta)
    GRABIX_CHECK(meta_fasta)
    ch_versions = ch_versions.mix(GRABIX_CHECK.out.versions)

    input_fasta = params.input
    if ( (params.input.endsWith(".gz")) && (GRABIX_CHECK.out.compress_bgzip == 'no')) {
        // Fasta is gzipped, but not bgzipped
        // This will cause WFMASH to complain
        // We will decompress it here, this way it is treated as a raw Fasta
        TABIX_BGZIP(meta_fasta)

        // Update the input Fasta with the new path
        input_fasta = TABIX_BGZIP.out.output
    }

    if (input_fasta.endsWith(".gz")) {
        meta_fasta = tuple([ id:fasta_file_name ], fasta)

        // For now we assume it was bgzip. If not WFMASH will complain instantly anyhow.
        if (!fai_path.exists() || !gzi_path.exists()) { // the assumption is that none of these files exist if only one does not exist
            SAMTOOLS_FAIDX(meta_fasta, [[],[]])
            fai = SAMTOOLS_FAIDX.out.fai
            gzi = SAMTOOLS_FAIDX.out.gzi
            ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        } else {
            fai = Channel.of([ [ id:fasta_file_name ], fai_path ])
            gzi = Channel.of([ [ id:fasta_file_name ], gzi_path ])
        }
        ch_fasta = meta_fasta
    } else {
        if (input_fasta.endsWith("fa")) {
            fasta_file_name = fasta_file_name.substring(0, fasta_file_name.length() - 3)
        } else {
            if (input_fasta.endswith("fasta")) {
                fasta_file_name = fasta_file_name.substring(0, fasta_file_name.length() - 6)
            } else { // we assume "fna" here
                fasta_file_name = fasta_file_name.substring(0, fasta_file_name.length() - 4)
            }
        }
        meta_fasta = tuple([ id:fasta_file_name ], fasta)
        TABIX_BGZIP(meta_fasta)
        ch_fasta = TABIX_BGZIP.out.output
        SAMTOOLS_FAIDX(ch_fasta, [[],[]])
        gzi = SAMTOOLS_FAIDX.out.gzi
        fai = SAMTOOLS_FAIDX.out.fai
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions)
    }

    emit:
    fasta = ch_fasta         // channel: [ val(meta), [ fasta ] ]
    fai = fai                // channel: [ val(meta), fasta.fai ]
    gzi = gzi                // channel: [ val(meta), fasta.gzi ]
    versions = ch_versions   // channel: [ versions.yml ]
}
