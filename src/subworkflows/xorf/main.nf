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
include { CHUNKER as UNMASKED_CHUNKER } from '../../modules/chunker/main.nf'
include { CONCAT }       from '../../modules/concat/main.nf'
include { CONCAT as CONCAT_RAW }   from '../../modules/concat/main.nf'
include { CONCAT as CONCAT_RENAMED }   from '../../modules/concat/main.nf'
include { CONCAT as CONCAT_RENAMED_RAW }   from '../../modules/concat/main.nf'
include { GENOMEMASK_SELENO } from '../../modules/genomemask/seleno/main.nf'
include { GENEPRED_LINT } from '../../modules/genepred/lint/main.nf'
include { DETACH_DUPLICATES } from '../../modules/detach/main.nf'
include { ISOTOOLS_TRUNCATION_DETECTOR } from '../../modules/isotools/utr/main.nf'
include { STRIP_OCCURRENCES as STRIP_TRUNCATIONS } from '../../modules/strip/main.nf'
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_UNMASKED } from '../../modules/bedtools/intersect/main.nf'
include { RENAME_PREDICTIONS as RENAME_PREDICTIONS_RAW } from '../../modules/rename/main.nf'
include { RENAME_PREDICTIONS } from '../../modules/rename/main.nf'

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
      output_dir     // path
      chunk_size     // int
      samba_weights  // channel: [ meta, path ]
      predict_keep_raw     // boolean
      selenocysteine_sites // path
      skip_netstart        // boolean
      rename_deactivate    // boolean

    main:
      def ch_regions  = regions
      def ch_sequence = sequence
      def ch_database = database
      def chunkSize   = chunk_size ?: 20

      def localHash = randomHash()

      def ch_versions = Channel.empty()

      // Linting //////////////////////////////////////////////////////////

      GENEPRED_LINT(
        ch_regions
      )

      // Chunking /////////////////////////////////////////////////////////

      ch_unmasked_chunked_regions = Channel.empty()
      ch_unmasked_chunked_sequences = Channel.empty()
      if (selenocysteine_sites) {
          ch_selenocysteine_sites = Channel.fromPath(selenocysteine_sites, checkIfExists: true).map { it -> [ [id: it.baseName], it ] }

          GENOMEMASK_SELENO(
              ch_sequence.map { it -> [ [id: it.baseName], it ] },
              ch_selenocysteine_sites
          )

          CHUNKER(
              ch_regions.map { meta, file -> [ meta + [ masked: true, hash: localHash ], file ] },
              GENOMEMASK_SELENO.out.fasta.first(),
              chunkSize,
          )

          BEDTOOLS_INTERSECT_UNMASKED(
              ch_regions,
              ch_selenocysteine_sites
          )
          UNMASKED_CHUNKER(
              BEDTOOLS_INTERSECT_UNMASKED.out.bed
                .map { meta, file -> [ meta + [ chr:'UNMSK', masked: false, hash: localHash ], file ] },
              ch_sequence.map { it -> [ [ id: it.baseName ], it ] },
              chunkSize,
          )
          ch_unmasked_chunked_regions = UNMASKED_CHUNKER.out.chunked_regions
          ch_unmasked_chunked_sequences = UNMASKED_CHUNKER.out.chunked_sequences

          ch_versions = ch_versions.mix(GENOMEMASK_SELENO.out.versions)
          ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT_UNMASKED.out.versions)
      } else {
          CHUNKER(
              ch_regions,
              ch_sequence.map { it -> [ [ id: it.baseName, hash: localHash ], it ] },
              chunkSize,
          )
      }

      // Joining /////////////////////////////////////////////////////////

      CHUNKER.out.chunked_regions
          .mix(ch_unmasked_chunked_regions)
          .flatMap { meta, region -> 
              def regions = region instanceof List ? region : [region]
              regions.collect { it ->
                [ meta + [ name: it.baseName ], it] }
              }
          .join(
              CHUNKER.out.chunked_sequences
                .mix(ch_unmasked_chunked_sequences)
                .flatMap { meta, fasta -> 
                    def fas = fasta instanceof List ? fasta : [fasta]
                    fas.collect { it ->
                      [ meta + [ name: it.baseName ], it] }
                }
          )
          .set { ch_pairs }

      // Prediction ///////////////////////////////////////////////////////

      GET_CANDIDATES(
          ch_pairs,
          ch_database,
          samba_weights,
          skip_netstart
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
              return tuple([ id: groupKey, name: metas[0].id, chr: metas[0].chr, masked: metas[0].masked, hash: localHash ], beds, tsvs)
          }
          .set { ch_all }

      // Renaming /////////////////////////////////////////////////////////

      ch_raw_renamed = Channel.empty()
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
                  return tuple([ name: metas[0].id, chr: metas[0].chr, id: groupKey, masked: metas[0].masked, hash: localHash ], beds, tsvs)
              }
              .set { ch_raw }

          CONCAT_RAW(
              ch_raw
          )

          if (!rename_deactivate) {
              RENAME_PREDICTIONS_RAW(
                  CONCAT_RAW.out.bed,
                  CONCAT_RAW.out.tsv
              )

             RENAME_PREDICTIONS_RAW.out.files
              .map { metas, beds, tsvs ->
                  [beds, tsvs, metas]
              }
              .collect(flat: false)
              .map { pairs ->
                  def bedFiles = pairs.collectMany { pair ->
                      pair[0] instanceof Collection ? pair[0] : [pair[0]]
                  }

                  def tsvFiles = pairs.collectMany { pair ->
                      pair[1] instanceof Collection ? pair[1] : [pair[1]]
                  }

                  def metadataFiles = pairs.collectMany { pair ->
                      pair[2] instanceof Collection ? pair[2] : [pair[2]]
                  }

                  def newId = metadataFiles[0].id 
                  return [
                      [ id: newId, name: localHash ],
                      bedFiles,
                      tsvFiles
                  ]
              }
              .set { ch_raw_renamed }
          }
      }
      
      // Concatenation /////////////////////////////////////////////////////

      CONCAT(
          ch_all
      )

      // Renaming /////////////////////////////////////////////////////////

      ch_full_length_with_duplicates = Channel.empty()
      if (!rename_deactivate) {
          RENAME_PREDICTIONS(
              CONCAT.out.bed,
              CONCAT.out.tsv
          )

          RENAME_PREDICTIONS.out.files
              .map { metas, beds, tsvs ->
                  [beds, tsvs, metas]
              }
              .collect(flat: false)
              .map { pairs ->
                  def bedFiles = pairs.collectMany { pair ->
                      pair[0] instanceof Collection ? pair[0] : [pair[0]]
                  }

                  def tsvFiles = pairs.collectMany { pair ->
                      pair[1] instanceof Collection ? pair[1] : [pair[1]]
                  }

                  def metadataFiles = pairs.collectMany { pair ->
                      pair[2] instanceof Collection ? pair[2] : [pair[2]]
                  }

                  def newId = metadataFiles[0].id
                  return [
                      [ id: newId, name: localHash ],
                      bedFiles,
                      tsvFiles
                  ]
              }
              .set { ch_renamed }

          // INFO: process uses 
          // > [meta.id, meta.name].findAll { it }.join('.') as new prefix
          CONCAT_RENAMED(
              ch_renamed
          )
          ch_full_length_with_duplicates = CONCAT_RENAMED.out.files.map { meta, bed, tsv -> tuple(meta, bed) }

          CONCAT_RENAMED_RAW(
              ch_raw_renamed
          )
      } else {
          ch_full_length_with_duplicates = CONCAT.out.files.map { meta, bed, tsv -> tuple(meta, bed) }
      }

      // Duplicates ///////////////////////////////////////////////////////

      DETACH_DUPLICATES( ch_full_length_with_duplicates )
      ch_full_length_transcripts = DETACH_DUPLICATES.out.pass
      ch_duplicates = DETACH_DUPLICATES.out.duplicates

      // Truncation ///////////////////////////////////////////////////////

      ISOTOOLS_TRUNCATION_DETECTOR(
          ch_full_length_transcripts
      )
      STRIP_TRUNCATIONS(
          ch_full_length_transcripts,
          ISOTOOLS_TRUNCATION_DETECTOR.out.descriptor,
          "TRUNCATED"
      )

      // Counts ///////////////////////////////////////////////////////////

      PREDICT_ORFS.out.counts
      .map { meta, initial, transaid, ns_td, tai, blast, samba, all, unique, kept ->
          def line = "${meta.id}@${meta.chr}\t${initial}\t${transaid}\t${ns_td}\t${tai}\t${blast}\t${samba}\t${all}\t${unique}\t${kept}"
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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def randomHash() {
    def chars = ('0'..'9') + ('a'..'z') + ('A'..'Z')
    def random = new Random()
    return (1..6).collect { chars[random.nextInt(chars.size())] }.join()
}
