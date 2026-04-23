// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GET_CANDIDATES — Runs translation, RNASamba, NetStart, TransAID, BLAST to get ORF candidates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { TRANSLATION } from '../../modules/translation/main.nf'
include { RNASAMBA }    from '../../modules/rnasamba/main.nf'
include { NETSTART }    from '../../modules/netstart/main.nf'
include { TRANSAID }    from '../../modules/transaid/main.nf'
include { BLAST }       from '../../modules/blast/main.nf'
include { JOIN as JOIN_NETS }   from '../../modules/join/main.nf'
include { WGET as WGET_SAMBA_WEIGHTS } from '../../modules/wget/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GET_CANDIDATES {
    take:
      ch_pairs
      ch_database
      samba_weights

    main:
      ch_versions = Channel.empty()

      TRANSLATION(ch_pairs)

      RNASAMBA(
        ch_pairs,
        samba_weights
      )

      NETSTART(
        RNASAMBA.out.fasta,
        RNASAMBA.out.bed
      )

      TRANSAID(
        RNASAMBA.out.fasta,
        RNASAMBA.out.bed
      )

       NETSTART.out.netstart
        .join(TRANSAID.out.transaid)
        .join(RNASAMBA.out.bed)
        .set { ch_nets }

      JOIN_NETS(
          ch_nets
      )

      TRANSLATION.out.predictions
      .join(JOIN_NETS.out.net)
      .set { ch_pre_candidates }

      BLAST(
          ch_pre_candidates,
          ch_database
      )
      .blast
      .join(RNASAMBA.out.samba)
      .set { ch_candidates }

      BLAST.out.counts
      .join(NETSTART.out.count)
      .join(TRANSAID.out.count)
      .join(JOIN_NETS.out.count)
      .set { ch_counts }

      ch_versions = ch_versions.mix(TRANSLATION.out.versions)
      ch_versions = ch_versions.mix(RNASAMBA.out.versions)
      ch_versions = ch_versions.mix(BLAST.out.versions)

    emit:
      candidates = ch_candidates
      counts = ch_counts
      versions = ch_versions
}
