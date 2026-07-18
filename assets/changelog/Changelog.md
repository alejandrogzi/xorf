# Changelog

All notable changes to this project are documented below.

---

## [0.0.37] - 2026-07-18

### Breaking Changes
- **BLAST/PSAURON/NETSTART2 re-impl**: Refactored the BLAST subworkflow to integrate PSAURON alongside DIAMOND for ORF scoring. The `orf blast` command now runs both tools sequentially, appending PSAURON confidence scores and reference sequence identifiers to each BLAST record. The output format has been extended accordingly.
- **`--skip_netstart` flag**: Added the ability to bypass NETSTART2 execution entirely. When `--skip_netstart` is set, the pipeline skips the NETSTART2 step and fills the merged net file with zeroed values for netstart-derived columns. This addresses cases where NETSTART2 is unavailable or undesirable.
- **NETSTART2 made optional in `orf net`**: The `--netstart` argument is no longer required. When omitted, the net merger receives an empty netstart map, yielding the same zeroed-out behaviour.

### orf v0.0.25
- **PSAURON integration**: The `orf blast` subcommand now runs PSAURON (`psauron`) after DIAMOND. Each ORF prediction is annotated with an in-frame score from PSAURON, and the output includes a `psauron_score` column.
- **Diamond `sseqid` in output**: The DIAMOND output format now includes the `sseqid` (subject sequence ID) column, providing the reference protein identifier for each match.
- **`--keep-temp` flag**: Added the `--prefix` and `--keep-temp` flags to `orf blast`. The `--prefix` flag controls output file naming, while `--keep-temp` retains intermediate files (Diamond raw output, PSAURON output, PEP files) instead of cleaning them up.
- **Graceful empty BLAST handling**: The pipeline no longer panics when DIAMOND produces zero predictions. An empty BLAST table is handled gracefully, preventing hard crashes on sequences with no database hits.
- **Container update**: The BLAST container now pins `diamond=2.2.4` and includes `psauron=1.1.3`.

### Infrastructure
- The `test` profile now uses `custom_database` instead of `database` for consistency.
- The `database` parameter in the `test` profile no longer points to the test dataset directly; instead, `custom_database` is used.
- Fixed a bug where `params.custom_database` was not being resolved correctly in the `GUNZIP_DATABASE` and direct `fromPath` branches.

## [0.0.36] - 2026-07-18

### Breaking Changes
- **CDS mapping fix for edge start/end cases**: Fixed the `map_cds` utility in `orf` (v0.0.24) to correctly handle edge cases where the CDS start or end falls at the boundaries of the alignment. This resolves incorrect ORF classifications for transcripts with start or stop codons exactly at exon junctions.
- **Default protein database retrieval**: Implemented automatic retrieval of the default UniProt/SwissProt protein database when no custom database is provided. The pipeline downloads the database from Zenodo if absent, and caches it for subsequent runs.
- **`custom_database` port**: Exposed the `custom_database` parameter as a first-class input, allowing users to supply a local DIAMOND database instead of relying on the default download.
- **UNMASKED default branch with selenocysteines**: When selenocysteine sites are provided, the pipeline now routes through an UNMASKED default branch, ensuring that selenocysteine masking produces an output with a clean `unmasked` identifier.

### orf v0.0.24
- Inclusive `map_cds` logic updated to correctly process edge positions.

## [0.0.35] - 2026-06-07

### Breaking Changes
- **Inclusive `map_cds`**: Updated `orf` v0.0.23 to make the CDS mapping function inclusive, fixing off-by-one errors that caused certain ORF boundaries to be misclassified.

### Chores
- BLAST processes are now submitted as array jobs for better scalability.

## [0.0.34] - 2026-06-05

### Breaking Changes
- **Net bug fix for non-netstart Transaid predictions**: Fixed a bug in `orf net` v0.0.22 where the net merger would produce incorrect output for transcripts lacking NETSTART2 predictions. Transaid-only predictions now merge correctly into the final net file.

### Features
- Added additional configuration options.
- Fixed documentation typos.
- Selenocysteine masking now outputs FASTA files directly.

## [0.0.33] - 2026-05-26

### Breaking Changes
- **xorf re-implementation [pre-release]**: Major internal rework of the xorf pipeline to improve modularity and data flow.

### Chores
- Updated README documentation.
- Left clear watermarks per job for traceability.

## [0.0.32] - 2026-05-22

### Breaking Changes
- **Detach/Split/Truncation + Arrays**: Implemented support for region detaching, splitting, and truncation. Added array job support for BLAST and other compute-intensive steps, enabling parallel execution across genomic regions.

### Fixes
- Updated container registry paths.
- Fixed version printing.

## [0.0.31] - 2026-05-09

### Breaking Changes
- **Genepred lint + docs**: Added gene prediction linting utilities and expanded documentation.
- **GenomeMask + Selenocysteine masking**: Implemented `genomemask` for selenocysteine masking.
- **Empty output handling**: Empty outputs from pipeline processes are now properly ignored instead of causing downstream failures.
- **Chunker input as tuple channel**: The chunker subworkflow now accepts sequence input as a tuple channel for better metadata propagation.

### Features
- Added `stageIn`/`stageOut` mode configuration options.

### Fixes
- Fixed ghost output emission from the predict step.
- Fixed selenocysteine channel consumption to avoid premature channel closure.
- Fixed WGET call for emitted channel; added support for local database paths.

## [0.0.30] - 2026-04-28

### Breaking Changes
- **Iso-chunk v0.0.21**: Fixed FASTA input handling. Added documentation and a `.gz` handler to prevent future bugs with compressed inputs.
- Switched to `flate2` default backend (`miniz_oxide`) for improved compatibility.

## [0.0.29] - 2026-04-23

### Breaking Changes
- **Chunker `--prefix` implementation**: Added `--prefix` parameter to `orf chunk` (v0.0.20), allowing configurable output file naming.

## [0.0.28] - 2026-04-21

### Breaking Changes
- **Chunker v0.0.19**: Fixed byte-seeking logic to correctly delete empty chunks. Added process tags with names for better traceability.

### Features
- Added `predict_keep_raw` extra argument to preserve raw prediction outputs.

## [0.0.27] - 2026-04-10

### Breaking Changes
- **Predict v0.0.17**: Implemented `concat_raw` for raw output concatenation. Added prediction ID support. Raw TSV and BED output available under the `--keep-raw` (`-K`) flag.

## [0.0.26] - 2026-04-02

### Breaking Changes
- **Chunker v0.0.18**: Implemented graceful handling of empty files when `--ignore-errors` is active. Fixed underflowing error when both flanks are applied to small sequences.

## [0.0.25] - 2026-03-25

### Breaking Changes
- **Predict v0.0.15**: New fixed XGBoost model for TranslationAI predictions, improving accuracy.
- **NetStart offline tokenizer**: Replaced the online tokenizer with an offline version, avoiding maximum request limit errors when queue sizes exceed 500.

## [0.0.24] - 2026-03-19

- Updated prediction model.
- Bug fixes.

## [0.0.23] - 2026-03-12

- Fixed predict headers to include net columns.
- Bug fixes and stability improvements.

## [0.0.22] - 2026-03-05

- Fixed typo in netstart header within predictions.
- Predict v0.0.10.

## [0.0.21] - 2026-02-26

- Fixed ORF coordinates across the board.
- Net shared ORF line tab spacing fix.

## [0.0.20] - 2026-02-19

- Fixed net Dockerfile and NetStart model directory.
- Stability improvements.

## [0.0.19] - 2026-02-12

- Fixed DIAMOND ghost records.
- Predict v0.0.8.

## [0.0.18] - 2026-02-05

- Fixed DIAMOND ghost records (re-fix).
- Predict v0.0.7.

## [0.0.17] - 2026-01-29

- Predict v0.0.6: fixed TranslationAI prefix splitting.
- Bug fixes.

## [0.0.16] - 2026-01-22

- **ORF blast v0.0.13**: Fixed unique TranslationAI NMD type assignment.

## [0.0.15] - 2026-01-15

- NMD detector fixed for negative strand coordinates.
- ORF blast v0.0.12.

## [0.0.14] - 2026-01-08

- Start/stop codon comparison made case-insensitive (ATCGatcgNn handling).

## [0.0.13] - 2026-01-01

- ATCGatcgNn handling across all sequence processing steps.

## [0.0.12] - 2025-12-25

- Detached candidate generation from prediction step.

## [0.0.11] - 2025-12-18

- Invalid codons accepted as X for downstream tolerance.

## [0.0.10] - 2025-12-11

- Fixed chunking slicing logic for edge cases.

## [0.0.9] - 2025-12-04

- Verbose error reporting for TranslationAI failures.

## [0.0.8] - 2025-11-27

- RNAsamba catching empty sequences gracefully.

## [0.0.7] - 2025-11-20

- Fixed exon slicing logic in the chunker.

## [0.0.6] - 2025-11-13

- Fixed chunking bug on extending exonic end.
- ORF v0.0.2.

## [0.0.5] - 2025-11-06

- Sequence slicing logic changed; added `--ignore-errors` flag.

## [0.0.4] - 2025-10-30

- Container builder CI setup.
- Initial release infrastructure.
