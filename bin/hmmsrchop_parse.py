#!/usr/bin/env python

"""
Parse the output of hmmsearch hmmsearch.out file and convert it to JSON format.

Usage:
    hmmsrchop.py parse <hmmsearch.out> [options]
    cat <hmmsearch.out> | hmmsrchop.py parse [options]

Required Arguments:
    <hmmsearch.out>   Path to the hmmsearch hmmsearch.out file.

Options:
    -h, --help     Show this help message and exit.
"""

import argparse
from dataclasses import (
    asdict,
    dataclass,
    field,
)
import json
import os
import pickle
from typing import (
    Any,
    Dict,
    Generator,
    List,
)
import sys

from Bio.SearchIO import parse
from Bio.SearchIO._model import QueryResult

from hmmsrchop_utils import (
    eprint,
)


SCHEMA_VERSION = "0.1.0"


@dataclass
class AlignmentRepresentation:
    profile_consensus: str
    match_annotation: str
    sequence: str
    posterior_prob: str

    def to_json(self) -> Dict[str, Any]:
        """
        Convert the AlignmentRepresentation object to a JSON-compatible dictionary.
        """
        return asdict(self)

    @classmethod
    def from_json(cls, json_data: Dict[str, Any]) -> "AlignmentRepresentation":
        """
        Create an AlignmentRepresentation object from a JSON-compatible dictionary.
        """
        return cls(**json_data)


@dataclass
class Alignment:
    alignment_id: str
    independent_evalue: float
    conditional_evalue: float
    bitscore: float
    profile_start: int
    profile_end: int
    envelope_start: int
    envelope_end: int
    sequence_start: int
    sequence_end: int
    alignment_len: int
    alignment_repr: AlignmentRepresentation
    significant: bool = False
    best_evalue: bool = False
    best_bitscore: bool = False
    best_aln_len: bool = False

    def to_json(self) -> Dict[str, Any]:
        """
        Convert the Alignment object to a JSON-compatible dictionary.
        """
        d = asdict(self)
        d["alignment_repr"] = self.alignment_repr.to_json()
        return d

    @classmethod
    def from_json(cls, d: dict) -> "Alignment":
        """
        Create an Alignment object from a JSON-compatible dictionary.
        """
        d["alignment_repr"] = AlignmentRepresentation.from_json(d["alignment_repr"])
        return cls(**d)


@dataclass
class ProfileHit:
    name: str
    accession: str
    description: str
    evalue: float
    bitscore: float
    profile_len: int
    alignments: List[Alignment] = field(default_factory=list)
    significant: bool = False
    best_evalue: bool = False
    best_bitscore: bool = False

    def flag_best_alignment(self):
        """
        Flag the alignment with the lowest E-value, the highest bit-score
        and the longest alignment length.

        Raises:
            IndexError: If there are no alignments for the hit.
        """
        if not self.alignments:
            raise IndexError(f"No alignments found for hit {self.accession}.")

        # E-value
        best_eval_index = min(
            range(len(self.alignments)),
            key=lambda i: self.alignments[i].conditional_evalue,
        )
        self.alignments[best_eval_index].best_evalue = True

        # Bit score
        best_bitscore_index = max(
            range(len(self.alignments)), key=lambda i: self.alignments[i].bitscore
        )
        self.alignments[best_bitscore_index].best_bitscore = True

        # Alignment length
        best_aln_len_index = max(
            range(len(self.alignments)), key=lambda i: self.alignments[i].alignment_len
        )
        self.alignments[best_aln_len_index].best_aln_len = True

    def to_json(self) -> Dict[str, Any]:
        """
        Convert the Hit object to a JSON-compatible dictionary.
        """
        d = asdict(self)
        d["alignments"] = [alignment.to_json() for alignment in self.alignments]
        return d

    @classmethod
    def from_json(cls, d: dict) -> "ProfileHit":
        """
        Create a Hit object from a JSON-compatible dictionary.
        """
        d["alignments"] = [
            Alignment.from_json(alignment) for alignment in d.get("alignments", [])
        ]
        return cls(**d)


@dataclass
class TargetSequence:
    """
    Stores the results associated with a single target sequence
    """

    sequence_id: str
    # sequence_len: int
    hits: List[ProfileHit] = field(default_factory=list)
    significant: bool = False
    schema_version: str = SCHEMA_VERSION

    def flag_best_hit(self):
        """
        Flag the hit with the lowest E-value and the one with the highest bit-score.

        Raises:
            IndexError: If there are no hits for the sequence.
        """
        if not self.hits:
            raise IndexError(f"No hits found for sequence {self.sequence_id}.")

        # E-value
        best_eval_index = min(range(len(self.hits)), key=lambda i: self.hits[i].evalue)
        self.hits[best_eval_index].best_evalue = True

        # Bit score
        best_bitscore_index = max(
            range(len(self.hits)), key=lambda i: self.hits[i].bitscore
        )
        self.hits[best_bitscore_index].best_bitscore = True

    def to_json_str(self) -> str:
        """
        Convert the TargetSequence object to a JSON-compatible dictionary.
        """
        d = asdict(self)
        d["hits"] = [hit.to_json() for hit in self.hits]
        return json.dumps(d, default=str)

    @classmethod
    def from_json_str(cls, json_data: str) -> "TargetSequence":
        """
        Create a TargetSequence object from a JSON-compatible dictionary.
        """
        d = json.loads(json_data)
        hits_data = d.pop("hits", [])
        hits = [ProfileHit.from_json(hit_data) for hit_data in hits_data]
        return cls(hits=hits, **d)

    def to_pickle_bytes(self) -> bytes:
        """
        Convert the TargetSequence object to a pickle byte stream.
        """
        return pickle.dumps(self)

    @classmethod
    def from_pickle_bytes(cls, pickle_data: bytes) -> "TargetSequence":
        """
        Create a TargetSequence object from a pickle byte stream.
        """
        return pickle.loads(pickle_data)


def id_generator(prefix: str = "DA_") -> Generator[str, None, None]:
    counter = 1
    while True:
        yield f"{prefix}{str(counter)}"
        counter += 1


def setup_argparse() -> argparse.ArgumentParser:
    """
        Sets up the argparse instance for command-line arguments.
    Returns:
        argparse.ArgumentParser: Configured ArgumentParser instance.
    """
    parser = argparse.ArgumentParser(
        add_help=False,
    )

    # Required arguments
    parser.add_argument(
        "hmmsearchout",
        metavar="<hmmsearch.out>",
        type=str,
        nargs="?",
    )

    parser.add_argument("-h", "--help", action="store_true", default=False)

    return parser


def setup_config(args: List[str]) -> argparse.Namespace:
    """
    Setup configuration for the script.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = setup_argparse()
    config = parser.parse_args(args)

    if len(sys.argv) == 1 and sys.stdin.isatty():
        print(__doc__)
        sys.exit(1)

    if config.help:
        print(__doc__)
        sys.exit(0)

    if not config.hmmsearchout:
        parser.print_help()
        eprint("Error: Missig required argument <hmmsearch.out")
        raise ValueError

    if not os.path.exists(config.hmmsearchout):
        eprint("Error: File not found: {}".format(config.hmmsearchout))
        raise FileNotFoundError

    return config


def parse_hmmsearch_out(hmmsearch_out: str) -> Generator:
    """
    Patse hmmsearch output file. Return a generator of query results.

    Args:
        hmmsearch_out (str): Path to hmmsearch.out

    Returns:
                   Generator of query results.
                   Each query result is a Bio.SearchIO object
                   For more information, see:
                     <https://biopython.org/docs/1.75/api/Bio.SearchIO.html>

    Raises:
        AssertionError: If the file does not exist.
    """

    if not os.path.isfile(hmmsearch_out):
        raise AssertionError("File does not exist: {}".format(hmmsearch_out))

    queryresult_generator = parse(hmmsearch_out, "hmmer3-text")

    return queryresult_generator


def parse_query_results(
    query_results: Generator,
    id_generator: Generator[str, None, None] = id_generator(),
) -> List[TargetSequence]:
    """
    Parse a single query result from the hmmsearch output.

    Args:
        query_result (QueryResult): A Bio.SearchIO QueryResult object.
        id_generator (Generator): A generator for unique IDs.

    Returns:
        TargetSequence: A TargetSequence object containing the parsed data.
    """
    sequences = {}
    alignment_id_generator = id_generator

    # The top level object is a QueryResult (from biopython)
    # In the hmmsearch output, it contains information mostly about the profile
    # But we're structuring the data in a sequence-centric way
    for qr in query_results:
        if not qr._items:
            continue  # This profile did not score any significant aligments

        # For the hmmsearch output, a profile scores hits against target sequences
        for hit in qr:
            sequence_id = hit.id

            # Do we have a new sequence?
            if sequence_id not in sequences:
                sequences[sequence_id] = TargetSequence(
                    sequence_id=sequence_id,
                )

            # Do we a new hit or is just a new alignment?
            hit_names = [h.name for h in sequences[sequence_id].hits]
            if qr.id not in hit_names:
                sequences[sequence_id].hits.append(
                    ProfileHit(
                        name=qr.id,
                        accession=""
                        if "accession" not in qr.__dict__
                        else qr.accession,
                        description=hit.description,
                        evalue=hit.evalue,
                        bitscore=hit.bitscore,
                        profile_len=qr.seq_len,
                        significant=hit.is_included,
                    )
                )

            for alignment in hit:
                # Is this possible?
                if len(alignment._items) > 1:
                    eprint(
                        "Error: More than one alignment found for hit {}.".format(
                            hit.id
                        )
                    )
                    eprint(
                        "Turns out we can have more than one alignment fragment for a single alignment."
                    )
                    sys.exit(1)

                alignment = Alignment(
                    alignment_id=next(alignment_id_generator),
                    independent_evalue=alignment.evalue,
                    conditional_evalue=alignment.evalue_cond,
                    bitscore=alignment.bitscore,
                    profile_start=alignment.env_start,
                    profile_end=alignment.env_end,
                    envelope_start=alignment.env_start,
                    envelope_end=alignment.env_end,
                    sequence_start=alignment._items[0].hit_start,
                    sequence_end=alignment._items[0].hit_end,
                    alignment_len=len(alignment._items[0].hit),
                    alignment_repr=AlignmentRepresentation(
                        profile_consensus=alignment._items[0].query.seq,
                        match_annotation=alignment._items[0].aln_annotation[
                            "similarity"
                        ],
                        sequence=alignment._items[0].hit.seq,
                        posterior_prob=alignment._items[0].aln_annotation["PP"],
                    ),
                    significant=hit.is_included,
                )

                for h in sequences[sequence_id].hits:
                    if h.name == qr.id:
                        h.alignments.append(alignment)
                        break

    for sequence in sequences.values():
        sequence.flag_best_hit()

        for hit in sequence.hits:
            # As you may find in the `hmmsearch` output:
            # 'No individual domains that satisfy reporting thresholds (although complete target did)'
            if hit.alignments:
                hit.flag_best_alignment()

    return list(sequences.values())


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
        queryresult_generator = parse_hmmsearch_out(config.hmmsearchout)
    except AssertionError as e:
        eprint(f"Error: {e}")
        sys.exit(1)

    try:
        parsed_query_results = parse_query_results(queryresult_generator)
    except ValueError as e:
        eprint("Error: {}".format(e))
        exit(1)

    for parsed_query_result in parsed_query_results:
        try:
            print(parsed_query_result.to_json_str())
        # <https://docs.python.org/3/library/signal.html#note-on-sigpipe>
        except BrokenPipeError:
            devnull = os.open(os.devnull, os.O_WRONLY)
            os.dup2(devnull, sys.stdout.fileno())
            sys.exit(1)
