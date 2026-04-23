#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    XORF — Main subworkflow: chunks input, gets candidates, predicts ORFs, concatenates results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CHUNKER }      from '../../modules/chunker/main.nf'
include { CONCAT }       from '../../modules/concat/main.nf'
include { CONCAT as CONCAT_RAW }   from '../../modules/concat/main.nf'
include { GENOMEMASK_SELENO } from '../../modules/genomemask/seleno/main.nf'
include { GENEPRED_LINT } from '../../modules/genepred/lint/main.nf'

include { PREDICT_ORFS } from '../predict_orfs/main.nf'
include { GET_CANDIDATES } from '../candidates/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow XORF {
    take:
      regions        // [ [id:id, chr:chr] , file ]
      sequence       // [ file ]
      database       // [ file ]
      output_dir
      chunk_size
      samba_weights  // channel: [ meta, path ]
      predict_keep_raw
      selenocysteine_sites

    main:
      def ch_regions  = regions
      def ch_sequence = sequence
      def ch_database = database
      def chunkSize   = chunk_size ?: 20

      def ch_versions = Channel.empty()

      GENEPRED_LINT(
        ch_regions
      )

      if (selenocysteine_sites) {
          GENOMEMASK_SELENO(
              ch_sequence.map { it -> [ [id: it.baseName], it ] },
              Channel.fromPath(selenocysteine_sites, checkIfExists: true)
              .map { it -> [ [id: it.baseName], it ] }
          )

          CHUNKER(
              ch_regions,
              GENOMEMASK_SELENO.out.twobit.first(),
              chunkSize,
          )

          ch_versions = ch_versions.mix(GENOMEMASK_SELENO.out.versions)
      } else {
          CHUNKER(
              ch_regions,
              ch_sequence.map { it -> [ [ id: it.baseName ], it ] },
              chunkSize,
          )
      }

      CHUNKER.out.chunked_regions
          .flatMap { meta, region -> 
              def regions = region instanceof List ? region : [region]
              regions.collect { it ->
                [ meta + [name: it.baseName], it] }
              }
          .join(
              CHUNKER.out.chunked_sequences
                .flatMap { meta, fasta -> 
                    def fas = fasta instanceof List ? fasta : [fasta]
                    fas.collect { it ->
                      [ meta + [name: it.baseName], it] }
                }
          )
          .set { ch_pairs }

      GET_CANDIDATES(
          ch_pairs,
          ch_database,
          samba_weights
      )

      PREDICT_ORFS(
          GET_CANDIDATES.out.candidates,
          GET_CANDIDATES.out.counts,
          predict_keep_raw
      )

      PREDICT_ORFS.out.orfs
          .filter { meta, bed, tsv -> bed.size() > 0 }
          .map { meta, bed, tsv -> 
              def groupKey = "${meta.id}@${meta.chr}"
              tuple(groupKey, meta, bed, tsv)
          }
          .groupTuple()
          .filter { groupKey, metas, beds, tsvs -> !beds.isEmpty() }
          .map { groupKey, metas, beds, tsvs ->
              return tuple([ id: groupKey, name: metas[0].id, chr: metas[0].chr ], beds, tsvs)
          }
          .set { ch_all }

      if (predict_keep_raw) {
          PREDICT_ORFS.out.raw
              .filter { meta, bed, tsv -> bed.size() > 0 }
              .map { meta, bed, tsv -> 
                  def groupKey = "${meta.id}@${meta.chr}"
                  tuple(groupKey, meta, bed, tsv)
              }
              .groupTuple()
              .filter { groupKey, metas, beds, tsvs -> !beds.isEmpty() }
              .map { groupKey, metas, beds, tsvs ->
                  return tuple([ name: metas[0].id, chr: metas[0].chr, id: groupKey ], beds, tsvs)
              }
              .set { ch_raw }

          CONCAT_RAW(
              ch_raw
          )
      }

      CONCAT(
          ch_all
      )

      PREDICT_ORFS.out.counts
      .map { meta, initial, netstart, transaid, ns_td, tai, blast, samba, all, unique, kept ->
          def line = "${meta.id}@${meta.chr}\t${initial}\t${netstart}\t${transaid}\t${ns_td}\t${tai}\t${blast}\t${samba}\t${all}\t${unique}\t${kept}"
          return line
      }
      .collectFile(
        name: 'counts.tsv',
        newLine: true,
        storeDir: "${output_dir}/XORF_PIPELINE_INFO/XORF_COUNTS"
      )
      .set { ch_counts }

      ch_versions = ch_versions.mix(CONCAT.out.versions)
      ch_pipeline_versions = ch_versions
          .collectFile(
              name:      "xorf.versions.yml",
              storeDir:  "${output_dir}/XORF_PIPELINE_INFO",
              sort:      true,
              keepHeader: false,
              newLine:   true
          )

      ch_versions = ch_versions.mix(CHUNKER.out.versions)
      ch_versions = ch_versions.mix(GET_CANDIDATES.out.versions)
      ch_versions = ch_versions.mix(PREDICT_ORFS.out.versions)

    emit:
      files = CONCAT.out.files
      counts = ch_counts
      versions = ch_pipeline_versions
}
