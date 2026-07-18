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
      ch_pairs       // channel: [ meta, path ]
      ch_database    // channel: [ meta, path ]
      samba_weights  // path
      skip_netstart  // boolean

    main:
      ch_versions = Channel.empty()

      TRANSLATION(ch_pairs)

      RNASAMBA(
        ch_pairs,
        samba_weights
      )

      TRANSAID(
        RNASAMBA.out.fasta,
        RNASAMBA.out.bed
      )

      if (!skip_netstart) {
        NETSTART(
          RNASAMBA.out.fasta,
          RNASAMBA.out.bed
        )

        ch_netstart = NETSTART.out.netstart
  
        TRANSAID.out.transaid
          .join(RNASAMBA.out.bed)
          .join(ch_netstart)
          .set { ch_nets }

        ch_nets.map { meta, transaid, bed, netstart -> [ meta, transaid, bed ] }.set { ch_transaid }
        ch_nets.map { meta, transaid, bed, netstart -> [ meta, netstart ] }.set { ch_netstart }

        JOIN_NETS( ch_transaid, ch_netstart )
        ch_versions = ch_versions.mix(NETSTART.out.versions)
      } else {
        TRANSAID.out.transaid
          .join(RNASAMBA.out.bed)
          .set { ch_nets }

        JOIN_NETS( ch_nets, Channel.value([[:], []]) )
      }

      TRANSLATION.out.predictions
      .join(JOIN_NETS.out.net)
      .set { ch_pre_candidates }

      // INFO: BLAST.out.blast is a channel of [meta, bed, tsv]
      // INFO: RNASAMBA.out.samba is a channel of [meta, tsv]
      BLAST(
          ch_pre_candidates,
          ch_database
      )
      .blast
      .join(RNASAMBA.out.samba)
      .set { ch_candidates }

      BLAST.out.counts
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
