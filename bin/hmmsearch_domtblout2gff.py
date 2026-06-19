#!/usr/bin/env python3
"""
Convert hmmsearch --domtblout output to GFF3.

Assumes the protein targets were produced by esl-translate from a six-frame
translation of a genome. The ORF coordinates and frame encoded in the
esl-translate FASTA header are used to project the protein hit back onto the
genomic nucleotide coordinates.

References:
  - HMMER User Guide (domtblout layout, p. 70): http://eddylab.org/software/hmmer/Userguide.pdf
  - esl-translate header format (esl-translate(1) man page):
        source=<seqname>  coords=<start>..<end>  length=<n_aa>  frame=<1..6>
        Frames 1-3 are the top strand, frames 4-6 are the bottom strand.
        On the bottom strand, start > end (start is the higher coordinate).
        Coords are 1-based, inclusive, and exclude the stop codon.
  - GFF3 spec: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
        GFF3 coordinates are 1-based, inclusive on both ends.

Coordinate math
---------------
ali_from / ali_to are 1-based protein positions in the ORF.

Top strand (frame 1..3):
    The ORF spans nucleotides [orf_start, orf_end] (orf_start < orf_end).
    Protein position p occupies nucleotides:
        [orf_start + 3*(p-1), orf_start + 3*p - 1]
    So the hit spans:
        nt_start = orf_start + 3*(ali_from - 1)
        nt_end   = orf_start + 3*ali_to - 1

Bottom strand (frame 4..6):
    esl-translate reports orf_start > orf_end; orf_start is the higher
    genomic coordinate (where translation begins). Protein position p
    occupies nucleotides closer to orf_end as p grows:
        [orf_start - 3*p + 1, orf_start - 3*(p-1)]
    So the hit spans:
        nt_start = orf_start - 3*ali_to + 1     (lower genomic coordinate)
        nt_end   = orf_start - 3*(ali_from - 1) (higher genomic coordinate)
    GFF3 requires start <= end, which this satisfies.

Usage
-----
    hmmsearch_domtblout2gff.py input.domtblout > output.gff3
    hmmsearch_domtblout2gff.py input.domtblout -o output.gff3
    cat input.domtblout | hmmsearch_domtblout2gff.py - > output.gff3
"""
from __future__ import annotations

import argparse
import datetime
import sys
from dataclasses import dataclass
from typing import Iterable, TextIO


# ---------------------------------------------------------------------------
# domtblout parsing
# ---------------------------------------------------------------------------

# Number of fixed, whitespace-delimited columns in domtblout before the
# free-form "description of target" field. The description can contain
# whitespace, so we split with maxsplit=N_FIXED_COLS to keep it intact.
N_FIXED_COLS = 22


@dataclass(frozen=True)
class DomtblRow:
    """A single domain hit row parsed from hmmsearch --domtblout."""
    target_name: str
    query_name: str
    query_accession: str
    query_length: int
    domain_i_evalue: float
    domain_score: float
    domain_bias: float
    hmm_from: int
    hmm_to: int
    ali_from: int
    ali_to: int
    env_from: int
    env_to: int
    accuracy: float
    description: str  # the trailing free-form field (esl-translate metadata)


def parse_domtblout_line(line: str) -> DomtblRow:
    """Parse one non-comment line of a hmmsearch --domtblout file."""
    # Keep the description intact by limiting split count.
    fields = line.rstrip("\n").split(maxsplit=N_FIXED_COLS)
    if len(fields) < N_FIXED_COLS + 1:
        raise ValueError(
            f"Expected at least {N_FIXED_COLS + 1} whitespace-delimited fields, "
            f"got {len(fields)}: {line!r}"
        )

    return DomtblRow(
        target_name     = fields[0],
        query_name      = fields[3],
        query_accession = fields[4],
        query_length    = int(fields[5]),
        domain_i_evalue = float(fields[12]),
        domain_score    = float(fields[13]),
        domain_bias     = float(fields[14]),
        hmm_from        = int(fields[15]),
        hmm_to          = int(fields[16]),
        ali_from        = int(fields[17]),
        ali_to          = int(fields[18]),
        env_from        = int(fields[19]),
        env_to          = int(fields[20]),
        accuracy        = float(fields[21]),
        description     = fields[22],
    )


# ---------------------------------------------------------------------------
# esl-translate header parsing
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class OrfInfo:
    """Source-genome coordinates and frame of an esl-translate ORF."""
    source: str      # name of the source DNA sequence (genome contig)
    orf_start: int   # 1-based; > orf_end on reverse strand
    orf_end: int
    orf_length: int  # length in amino acids
    frame: int       # 1..3 = top strand, 4..6 = bottom strand


def parse_esl_translate_description(description: str) -> OrfInfo:
    """Parse the esl-translate metadata fields from the domtblout description.

    Expected format (whitespace-separated key=value pairs, possibly followed
    by the original source description text we don't care about):
        source=<name>  coords=<start>..<end>  length=<n>  frame=<n>  [trailing text]
    """
    fields = {}
    for token in description.split():
        if "=" not in token:
            continue  # ignore trailing free-form description text
        key, _, value = token.partition("=")
        fields[key] = value

    for required in ("source", "coords", "length", "frame"):
        if required not in fields:
            raise ValueError(
                f"Missing '{required}=' field in esl-translate description: "
                f"{description!r}"
            )

    try:
        coord_start_s, coord_end_s = fields["coords"].split("..")
        orf_start = int(coord_start_s)
        orf_end   = int(coord_end_s)
    except ValueError as exc:
        raise ValueError(
            f"Could not parse coords field {fields['coords']!r}: {exc}"
        ) from exc

    return OrfInfo(
        source     = fields["source"],
        orf_start  = orf_start,
        orf_end    = orf_end,
        orf_length = int(fields["length"]),
        frame      = int(fields["frame"]),
    )


# ---------------------------------------------------------------------------
# Coordinate projection: protein hit -> genomic nucleotide interval
# ---------------------------------------------------------------------------

def project_to_genome(
    ali_from: int,
    ali_to: int,
    orf: OrfInfo,
) -> tuple[int, int, str]:
    """Project a protein hit (ali_from..ali_to, 1-based incl.) onto the genome.

    Returns (nt_start, nt_end, strand) with GFF3 conventions:
        - 1-based, inclusive on both ends
        - nt_start <= nt_end regardless of strand
        - strand in {'+', '-'}

    Raises ValueError on invalid frame or inconsistent coordinates.
    """
    if not (1 <= ali_from <= ali_to <= orf.orf_length):
        raise ValueError(
            f"ali_from..ali_to ({ali_from}..{ali_to}) out of range for ORF of "
            f"length {orf.orf_length}"
        )

    if 1 <= orf.frame <= 3:
        # Top strand. orf_start < orf_end (we assert below).
        if orf.orf_start > orf.orf_end:
            raise ValueError(
                f"Frame {orf.frame} is top-strand but orf_start ({orf.orf_start}) "
                f"> orf_end ({orf.orf_end})"
            )
        nt_start = orf.orf_start + 3 * (ali_from - 1)
        nt_end   = orf.orf_start + 3 * ali_to - 1
        strand   = "+"
    elif 4 <= orf.frame <= 6:
        # Bottom strand. esl-translate sets orf_start > orf_end.
        if orf.orf_start < orf.orf_end:
            raise ValueError(
                f"Frame {orf.frame} is bottom-strand but orf_start ({orf.orf_start}) "
                f"< orf_end ({orf.orf_end})"
            )
        nt_start = orf.orf_start - 3 * ali_to + 1
        nt_end   = orf.orf_start - 3 * (ali_from - 1)
        strand   = "-"
    else:
        raise ValueError(f"Invalid frame {orf.frame!r}; expected 1..6")

    # Sanity-check that the projected interval falls inside the ORF span.
    orf_lo = min(orf.orf_start, orf.orf_end)
    orf_hi = max(orf.orf_start, orf.orf_end)
    if not (orf_lo <= nt_start <= nt_end <= orf_hi):
        raise ValueError(
            f"Projected nt interval [{nt_start}, {nt_end}] falls outside ORF span "
            f"[{orf_lo}, {orf_hi}] (ali_from={ali_from}, ali_to={ali_to}, "
            f"frame={orf.frame})"
        )

    return nt_start, nt_end, strand


# ---------------------------------------------------------------------------
# GFF3 emission
# ---------------------------------------------------------------------------

# Per the GFF3 spec, these characters MUST be percent-encoded when they appear
# in attribute keys or values: tab, newline, carriage return, ;, =, &, %, and
# control chars. We keep this minimal because the values we emit are bounded
# (no user-supplied free text in the keys).
_GFF3_ATTR_ESCAPE = str.maketrans({
    "\t": "%09",
    "\n": "%0A",
    "\r": "%0D",
    ";":  "%3B",
    "=":  "%3D",
    "&":  "%26",
    ",":  "%2C",
})


def _escape_attr(value: str) -> str:
    return value.translate(_GFF3_ATTR_ESCAPE)


def build_attributes(row: DomtblRow, orf: OrfInfo) -> str:
    """Build the GFF3 attributes (column 9) for one hit."""
    # Order chosen so ID comes first (GFF3 convention) and HMM-related fields
    # are grouped together.
    pairs = [
        ("ID",              row.target_name),
        ("Query",           row.query_name),
        ("QueryAccession",  row.query_accession),
        ("QueryLength",     str(row.query_length)),
        ("iEvalue",         f"{row.domain_i_evalue:g}"),
        ("DomainScore",     f"{row.domain_score:g}"),
        ("DomainBias",      f"{row.domain_bias:g}"),
        ("HMMFrom",         str(row.hmm_from)),
        ("HMMTo",           str(row.hmm_to)),
        ("AliFrom",         str(row.ali_from)),
        ("AliTo",           str(row.ali_to)),
        ("EnvFrom",         str(row.env_from)),
        ("EnvTo",           str(row.env_to)),
        ("Accuracy",        f"{row.accuracy:g}"),
        ("OrfStart",        str(orf.orf_start)),
        ("OrfEnd",          str(orf.orf_end)),
        ("OrfFrame",        str(orf.frame)),
        ("OrfLength",       str(orf.orf_length)),
    ]
    return ";".join(f"{k}={_escape_attr(v)}" for k, v in pairs)


def format_gff3_line(row: DomtblRow, orf: OrfInfo) -> str:
    """Build one GFF3 record from a domtbl row + its esl-translate metadata."""
    nt_start, nt_end, strand = project_to_genome(row.ali_from, row.ali_to, orf)
    return "\t".join((
        orf.source,                # seqid
        "hmmsearch",               # source
        "hmmsearch_hit",           # type
        str(nt_start),             # start
        str(nt_end),               # end
        f"{row.domain_score:g}",   # score
        strand,                    # strand
        ".",                       # phase (N/A for non-CDS features)
        build_attributes(row, orf),
    ))


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------

def _iso_timestamp() -> str:
    # Local time with tzinfo, ISO-8601 second precision.
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def convert(
    in_stream: TextIO,
    out_stream: TextIO,
    on_error: str = "fail",
) -> tuple[int, int]:
    """Stream-convert a domtblout file to GFF3.

    Comment/metadata lines from the *bottom* of the domtblout (the per-run
    summary block) are preserved verbatim at the end of the output, matching
    the original awk behaviour. The first three header lines (the column
    titles) are dropped.

    Args:
        in_stream:  open text stream with the domtblout content.
        out_stream: open text stream to write GFF3 to.
        on_error:   'fail' (default) raises on the first malformed line;
                    'skip' logs to stderr and continues.

    Returns:
        (n_hits_written, n_lines_skipped)
    """
    if on_error not in ("fail", "skip"):
        raise ValueError(f"on_error must be 'fail' or 'skip', got {on_error!r}")

    print("##gff-version 3", file=out_stream)
    print("# Generated from hmmsearch --domtblout by hmmsearch_domtblout2gff.py",
          file=out_stream)
    print(f"# Timestamp: {_iso_timestamp()}", file=out_stream)

    trailing_comments: list[str] = []
    n_written = 0
    n_skipped = 0

    for line_no, raw_line in enumerate(in_stream, start=1):
        line = raw_line.rstrip("\n")
        if not line.strip():
            continue

        if line.startswith("#"):
            # Drop the three-line header block at the top; keep the trailing
            # per-run metadata block.
            if line_no <= 3:
                continue
            trailing_comments.append(line)
            continue

        try:
            row = parse_domtblout_line(line)
            orf = parse_esl_translate_description(row.description)
            print(format_gff3_line(row, orf), file=out_stream)
            n_written += 1
        except (ValueError, IndexError) as exc:
            msg = f"line {line_no}: {exc}"
            if on_error == "fail":
                raise ValueError(msg) from exc
            print(f"WARN: skipping {msg}", file=sys.stderr)
            n_skipped += 1

    # Preserve the trailing metadata block at the end of the file so it is
    # not lost (it documents the hmmsearch invocation, version, etc.).
    for comment in trailing_comments:
        print(comment, file=out_stream)

    return n_written, n_skipped


def _build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Convert hmmsearch --domtblout (with esl-translate target "
                    "headers) to GFF3 in genomic coordinates.",
    )
    p.add_argument(
        "input",
        help="Path to the hmmsearch --domtblout file, or '-' to read stdin.",
    )
    p.add_argument(
        "-o", "--output",
        default="-",
        help="Output GFF3 path, or '-' for stdout (default: stdout).",
    )
    p.add_argument(
        "--on-error",
        choices=("fail", "skip"),
        default="fail",
        help="Behaviour on malformed lines (default: fail).",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_argparser().parse_args(argv)

    in_stream: TextIO
    out_stream: TextIO
    if args.input == "-":
        in_stream = sys.stdin
    else:
        in_stream = open(args.input, "r", encoding="utf-8")

    if args.output == "-":
        out_stream = sys.stdout
    else:
        out_stream = open(args.output, "w", encoding="utf-8")

    try:
        n_written, n_skipped = convert(in_stream, out_stream, on_error=args.on_error)
    finally:
        if in_stream is not sys.stdin:
            in_stream.close()
        if out_stream is not sys.stdout:
            out_stream.close()

    print(
        f"Wrote {n_written} hits"
        + (f" ({n_skipped} lines skipped)" if n_skipped else ""),
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
