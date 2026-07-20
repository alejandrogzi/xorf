#!/usr/bin/env python3
"""Rename ORF predictions in a BED12 file and its companion TSV records.

This script takes a BED12 file of predicted open reading frames (ORFs) and a
tab-separated (TSV) file holding additional metadata for each prediction, and
rewrites the shared ``id`` column of both files with human-readable,
information-rich names.

Names are built by combining an optional prefix, a root identifier, optional
tags, an optional ORF score, and an optional protein (BLAST) reference. The
resulting identifier has the general shape::

    [PREFIX#...]ROOT[#TAG]...[#SC{SCORE}][@PROTEIN]

By default the root identifier is the part of the original id before the ``@``
separator. With ``--rebase`` the entire original id is replaced by a short,
unique, content-derived hash so that identical predictions across runs receive
identical names.

Reserved characters ``#`` and ``@`` are not permitted inside user-supplied
prefixes or tags because they are used as the structural separators of the
generated names.

Example
-------
    rename_predictions.py \\
        --bed predictions.bed \\
        --tsv predictions.tsv \\
        --unmask --tag LIVER --append-orf-score \\
        --output-bed renamed.bed --output-tsv renamed.tsv

See :func:`parse_args` for the full list of command-line options.
"""

import argparse
import hashlib
import sys

import pandas as pd

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "alejandrxgzi@gmail.com"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.0.1"

BED_COLS = [
    "chrom",
    "chrom_start",
    "chrom_end",
    "id",
    "score",
    "strand",
    "thick_start",
    "thick_end",
    "item_rgb",
    "block_count",
    "block_sizes",
    "block_starts",
]
RESERVED = {"#", "@"}


def validate_token(value: str, label: str) -> str:
    """Reject tokens containing reserved separator characters.

    Prefixes and tags must not include ``#`` or ``@`` because those characters
    are used to separate the components of a generated name. ``label`` is used
    to produce a helpful error message identifying which argument was invalid.
    """
    if any(char in value for char in RESERVED):
        raise ValueError(f"{label} cannot contain '#' or '@': {value!r}")
    return value


def short_hashes(ids: list[str]) -> dict[str, str]:
    """Map each unique identifier to a short, collision-free hash.

    For every distinct ``id`` a 64-bit BLAKE2s digest is computed. The digest
    is then truncated to the smallest prefix length (between 4 and 16 hex
    characters) that keeps all generated hashes unique among the inputs, and
    prefixed with ``h`` to form the replacement root name.

    Raises
    ------
    ValueError
        If no prefix width between 4 and 16 yields a set of unique hashes,
        which can only happen with an impractically large number of inputs.
    """
    unique_ids = sorted(set(ids))
    digests = {
        value: hashlib.blake2s(value.encode(), digest_size=8).hexdigest()
        for value in unique_ids
    }

    for width in range(4, 17):
        values = [digest[:width] for digest in digests.values()]
        if len(values) == len(set(values)):
            return {key: f"h{digest[:width]}" for key, digest in digests.items()}

    raise ValueError("Could not generate unique short hashes")


def format_score(value: str) -> str:
    """Format an ORF score as a compact integer string.

    The value is parsed as a integer (0-100) and rendered with two significant digits
    (``%.2g``) so scores from the TSV are normalized consistently in the
    generated names.
    """
    try:
        return f"{float(value) * 100:.0f}"
    except ValueError as error:
        raise ValueError(f"Invalid ORF score: {value!r}") from error


def parse_args() -> argparse.Namespace:
    """Define and parse the command-line interface.

    Returns
    -------
    argparse.Namespace
        Parsed arguments. Available attributes:

        ``bed``
            Path to the input BED12 file.
        ``tsv``
            Path to the input TSV metadata file.
        ``output_bed``
            Destination path for the renamed BED12 file.
        ``output_tsv``
            Destination path for the renamed TSV file.
        ``deactivate``
            When set, copy the inputs unchanged to the outputs and exit.
        ``rebase``
            Replace the entire id with a short content hash instead of keeping
            the original root.
        ``append_orf_score``
            Append the ORF score (from ``score_column``) to each name.
        ``score_column``
            TSV column used as the ORF score source.
        ``append_protein_name``
            Append the BLAST reference id as the protein component.
        ``custom_prefix``
            Extra prefix prepended to every generated name.
        ``unmask``
            Prepend the ``UNMASK`` prefix to every generated name.
        ``tag``
            Zero or more tags appended to every generated name.
    """
    parser = argparse.ArgumentParser(
        description="Rename BED predictions and their TSV records."
    )
    parser.add_argument("--bed", required=True)
    parser.add_argument("--tsv", required=True)
    parser.add_argument("--output-bed", default="renamed.bed")
    parser.add_argument("--output-tsv", default="renamed.tsv")
    parser.add_argument("--deactivate", action="store_true")
    parser.add_argument("--rebase", action="store_true")
    parser.add_argument("--append-orf-score", action="store_true")
    parser.add_argument("--score-column", default="prob_coding")
    parser.add_argument("--append-protein-name", action="store_true")
    parser.add_argument("--custom-prefix")
    parser.add_argument("--unmask", action="store_true")
    parser.add_argument("--tag", action="append", default=[])
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    bed = pd.read_csv(args.bed, sep="\t", header=None, dtype=str, keep_default_na=False)
    if bed.shape[1] < len(BED_COLS):
        raise ValueError(f"Expected BED12, found {bed.shape[1]} columns")

    if bed.shape[1] > len(BED_COLS):
        print(f"WARN: BED12 has {bed.shape[1]} columns; dropping extra columns")
        bed.drop(columns=bed.columns[len(BED_COLS):], inplace=True)

    bed.columns = BED_COLS

    tsv = pd.read_csv(args.tsv, sep="\t", dtype=str, keep_default_na=False)
    if "id" not in tsv.columns:
        raise ValueError("TSV is missing the 'id' column")
    if tsv["id"].duplicated().any():
        raise ValueError("TSV contains duplicate IDs")

    bed_ids = set(bed["id"])
    tsv_ids = set(tsv["id"])
    if bed_ids != tsv_ids:
        missing_tsv = sorted(bed_ids - tsv_ids)[:5]
        missing_bed = sorted(tsv_ids - bed_ids)[:5]
        raise ValueError(
            f"BED/TSV IDs differ; missing from TSV: {missing_tsv}; "
            f"missing from BED: {missing_bed}"
        )

    if args.deactivate:
        bed.to_csv(args.output_bed, sep="\t", header=False, index=False)
        tsv.to_csv(args.output_tsv, sep="\t", index=False)
        return

    prefix_parts = []
    if args.unmask:
        prefix_parts.append("UNMASK")
    if args.custom_prefix:
        prefix_parts.append(validate_token(args.custom_prefix, "custom prefix"))

    tags = [validate_token(tag, "tag").upper() for tag in args.tag]
    hash_map = short_hashes(tsv["id"].tolist()) if args.rebase else {}
    records = tsv.set_index("id", drop=False).to_dict("index")

    if args.append_orf_score and args.score_column not in tsv.columns:
        raise ValueError(f"TSV is missing score column {args.score_column!r}")
    if args.append_protein_name and "blast_reference_id" not in tsv.columns:
        raise ValueError("TSV is missing the 'blast_reference_id' column")

    rename_map = {}
    for old_id, row in records.items():
        root, separator, existing_protein = old_id.partition("@")
        name = f"xORF#{hash_map[old_id]}" if args.rebase else root

        if prefix_parts:
            name = "#".join(prefix_parts + [name])
        if tags:
            name += "#" + "#".join(tags)
        if args.append_orf_score and row[args.score_column] != "":
            name += f"#SC{format_score(row[args.score_column])}"

        protein = row.get("blast_reference_id", "") if args.append_protein_name else ""
        if not args.rebase and not args.append_protein_name and separator:
            protein = existing_protein
        if protein:
            name += f"@{protein}"

        rename_map[old_id] = name

    bed["id"] = bed["id"].map(rename_map)
    tsv["id"] = tsv["id"].map(rename_map)

    bed.to_csv(args.output_bed, sep="\t", header=False, index=False)
    tsv.to_csv(args.output_tsv, sep="\t", index=False)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, pd.errors.ParserError) as error:
        sys.exit(f"ERROR: {error}")
