#!/usr/bin/env python

"""
Write the coordinates of the significant HMM profile alignments of each sequence to GFF format.

Usage:
    hmmsrchop.py togff <serialized HMMER3 output>
    hmmsrchop.py togff -h | --help
    cat <serialized HMMER3 output> | hmmsrchop.py togff

Required Arguments:
    <serialized HMMER3 output>  The serialized HMMER3 output file to convert to GFF format.

Modes of Operation:
    --mode <mode>               Conversion mode. Options are:
                                    all-alignments      All aligments from all profiles
                                    best-profile        All aligments of the best performing profile
                                    best-alignment      Single best alignment from the best performing profile
                                    best-per-profile    Single best alignment from each profile
                                [default: best-alignment]

Metrics for Selecting Best-Scoring Profiles:
    --metric <metric>           Metric to use for selecting the best scoring profile. Options are:
                                    bitscore            Bit score
                                    evalue              E-value
                                [default: evalue]

General Options:
    -h, --help                  Show this help message and exit.
        -v, --version               Show version information and exit.

Examples:
  # Convert best alignment per profile using E-value
  cat hmmsearch.out | hmmsrchop.py parse | hmmsrchop.py togff > best_alignments.gff

  # Convert all alignments from all profiles using E-value
  cat hmmsearch.out | hmmsrchop.py parse | hmmsrchop.py togff --mode all-alignments > all_alignments.gff

  # Single best alignment from best scoring profile by bitscore
  cat hmmsearch.out | hmmsrchop.py parse | hmmsrchop.py togff --mode best-profile --metric bitscore > best_profile.gff
"""

import argparse
from dataclasses import (
    asdict,
    dataclass,
    field,
)
from enum import Enum
import json
import os
import sys
from typing import (
    Dict,
    List,
    Tuple,
)

from hmmsrchop_parse import (
    Alignment,
    TargetSequence,
)
from hmmsrchop_utils import (
    eprint,
    parse_seqkit2translate_fasta_header,
    read_input,
)


class GffMode(Enum):
    ALL_ALIGNMENTS = "all-alignments"
    BEST_PROFILE = "best-profile"
    BEST_ALIGNMENT = "best-alignment"
    BEST_PER_PROFILE = "best-per-profile"


class GffMetric(Enum):
    BIT_SCORE = "bitscore"
    EVALUE = "evalue"


@dataclass
class GffFeature:
    """
    Represents a GFF feature with its attributes.

    Attributes:
        feature_id (str): Identifier of the feature. (Domain Alignment ID)
        seqid (str): Sequence identifier.
        source (str): Source of the feature.
        type_ (str): Type of the feature.
        start (int): Start position of the feature.
        end (int): End position of the feature.
        score (float | str): Score of the feature.
        strand (str): Strand of the feature.
        phase (int | str): Phase of the feature.
        attributes (Dict[str, str]): Additional attributes of the feature
    """

    feature_id: str
    seqid: str
    source: str
    type_: str
    start: int
    end: int
    score: float | str
    strand: str
    phase: int | str
    attributes: Dict[str, str]

    def __repr__(self) -> str:
        """
        Return a string representation of the GFF feature.

        Returns:
            str: String representation of the GFF feature.
        """
        return "\t".join(
            [
                self.seqid,
                self.source,
                self.type_,
                str(self.start),
                str(self.end),
                str(self.score),
                self.strand,
                str(self.phase),
                ";".join(f"{k}={v}" for k, v in self.attributes.items()),
            ]
        )

    @classmethod
    def from_profile_alignment(
        cls,
        profile_alignment: Alignment,
        sequence_id: str,
        profile_name: str,
    ) -> "GffFeature":
        """ """
        # the translation frame is not necessary as the translation start
        # and end positions are already in the correct frame
        parent_seq, tr_frame, tr_start, _ = parse_seqkit2translate_fasta_header(
            sequence_id
        )

        # WARNING: test the outcome if you touch the following lines
        alignment_start = tr_start + (profile_alignment.sequence_start * 3)
        alignment_end = (tr_start + (profile_alignment.sequence_end * 3)) - 1
        alignment_strand = "+" if tr_frame > 0 else "-"

        return GffFeature(
            feature_id=profile_alignment.alignment_id,
            seqid=parent_seq,
            source="hmmsrchop.py",
            type_="domain_alignment",
            start=alignment_start,
            end=alignment_end,
            score=profile_alignment.bitscore,
            strand=alignment_strand,
            phase=".",
            attributes={
                "featureID": profile_alignment.alignment_id,
                "profileName": profile_name,
                "conditionalEvalue": str(profile_alignment.conditional_evalue),
                "independentEvalue": str(profile_alignment.independent_evalue),
                "bitscore": str(profile_alignment.bitscore),
                "alignmentLength": str(profile_alignment.alignment_len),
            },
        )


def gff3_from_target_sequence(
    target_sequence: TargetSequence,
    mode: GffMode,
    metric: GffMetric,
) -> List[GffFeature]:
    """
    Convert the domain alignment results from a target sequence into GFF3 format.

    Args:
        target_sequence (TargetSequence): The target sequence containing the domain alignment results.
        mode (GffMode): The mode of operation for GFF conversion.
        metric (GffMetric): The metric used for selecting the best scoring profile.

    Returns:
        List[GffFeature]: A list of GFF features representing the domain alignments.
    """
    gff_features = []

    for profile_hit in target_sequence.hits:
        # Skip if no domain alignments
        if not profile_hit.alignments:
            continue
        if not profile_hit.significant:
            continue

        match mode:
            case GffMode.ALL_ALIGNMENTS:
                for alignment in profile_hit.alignments:
                    gff_features.append(
                        GffFeature.from_profile_alignment(
                            profile_alignment=alignment,
                            sequence_id=target_sequence.sequence_id,
                            profile_name=profile_hit.name,
                        )
                    )

            case GffMode.BEST_PROFILE:
                match metric:
                    case GffMetric.BIT_SCORE:
                        if profile_hit.best_bitscore:
                            for alignment in profile_hit.alignments:
                                gff_features.append(
                                    GffFeature.from_profile_alignment(
                                        profile_alignment=alignment,
                                        sequence_id=target_sequence.sequence_id,
                                        profile_name=profile_hit.name,
                                    )
                                )

                    case GffMetric.EVALUE:
                        if profile_hit.best_evalue:
                            for alignment in profile_hit.alignments:
                                gff_features.append(
                                    GffFeature.from_profile_alignment(
                                        profile_alignment=alignment,
                                        sequence_id=target_sequence.sequence_id,
                                        profile_name=profile_hit.name,
                                    )
                                )

            case GffMode.BEST_ALIGNMENT:
                match metric:
                    case GffMetric.BIT_SCORE:
                        if profile_hit.best_bitscore:
                            for alignment in profile_hit.alignments:
                                if alignment.best_bitscore:
                                    gff_features.append(
                                        GffFeature.from_profile_alignment(
                                            profile_alignment=alignment,
                                            sequence_id=target_sequence.sequence_id,
                                            profile_name=profile_hit.name,
                                        )
                                    )

                    case GffMetric.EVALUE:
                        if profile_hit.best_evalue:
                            for alignment in profile_hit.alignments:
                                if alignment.best_evalue:
                                    gff_features.append(
                                        GffFeature.from_profile_alignment(
                                            profile_alignment=alignment,
                                            sequence_id=target_sequence.sequence_id,
                                            profile_name=profile_hit.name,
                                        )
                                    )

            case GffMode.BEST_PER_PROFILE:
                match metric:
                    case GffMetric.BIT_SCORE:
                        for alignment in profile_hit.alignments:
                            if alignment.best_bitscore:
                                gff_features.append(
                                    GffFeature.from_profile_alignment(
                                        profile_alignment=alignment,
                                        sequence_id=target_sequence.sequence_id,
                                        profile_name=profile_hit.name,
                                    )
                                )

                    case GffMetric.EVALUE:
                        for alignment in profile_hit.alignments:
                            if alignment.best_evalue:
                                gff_features.append(
                                    GffFeature.from_profile_alignment(
                                        profile_alignment=alignment,
                                        sequence_id=target_sequence.sequence_id,
                                        profile_name=profile_hit.name,
                                    )
                                )

    return gff_features


def setup_argparse() -> argparse.ArgumentParser:
    """
    Sets up the argparse instance for command-line arguments.

    Returns:
    argparse.ArgumentParser: Configured ArgumentParser instance.
    """

    parser = argparse.ArgumentParser(
        add_help=False,
    )

    # Required Arguments
    parser.add_argument(
        "input",
        type=str,
        default="-",
        help="The serialized HMMER3 output file to convert to GFF format.",
    )

    # Modes of Operation
    parser.add_argument(
        "--mode",
        type=str,
        choices=[mode.value for mode in GffMode],
        default=GffMode.BEST_ALIGNMENT.value,
        help="Conversion mode. Options: {all-alignments, best-profile, best-alignment, best-per-profile}. [default: best-alignment]",
    )

    # Metrics for Selecting Best-Scoring Profiles
    parser.add_argument(
        "--metric",
        type=str,
        choices=[metric.value for metric in GffMetric],
        default=GffMetric.EVALUE.value,
        help="Metric to use for selecting the best scoring profile. Options: {bitscore, evalue}. [default: evalue]",
    )

    # General Options
    parser.add_argument("-h", "--help", action="store_true", default=False)

    return parser


def setup_config(args: List[str]) -> argparse.Namespace:
    """
    Setup configuration for the script.

    Returns:
        argparse.Namespace: Parsed command-line arguments.

    Raises:
        ValueError: If no input file is specified.
    """
    parser = setup_argparse()
    config = parser.parse_args(args)

    if len(sys.argv) == 1 and sys.stdin.isatty():
        print(__doc__)
        sys.exit(1)

    if config.help:
        print(__doc__)
        sys.exit(0)

    if not config.input:
        eprint("Error: No input file specified.")
        raise ValueError

    return config


def run(args: List[str]) -> None:
    try:
        config = setup_config(args)
    except (ValueError, FileNotFoundError):
        sys.exit(1)
    eprint(
        "Info: Running with the following configuration: "
        + ", ".join(f"{k}={v}" for k, v in config.__dict__.items())
    )

    try:
        file_handle = read_input(config.input)
    except FileNotFoundError:
        sys.exit(1)

    target_sequences = []
    for line in file_handle:
        target_sequence = TargetSequence.from_json_str(line)
        target_sequences.append(target_sequence)

    if not target_sequences:
        eprint("Error: No target sequences found.")
        sys.exit(1)

    for target_sequence in target_sequences:
        gff_features = gff3_from_target_sequence(
            target_sequence=target_sequence,
            mode=GffMode(config.mode),
            metric=GffMetric(config.metric),
        )

        for gff_feature in gff_features:
            try:
                print(gff_feature)
            except BrokenPipeError:
                devnull = os.open(os.devnull, os.O_WRONLY)
                os.dup2(devnull, sys.stdout.fileno())
                sys.exit(1)
