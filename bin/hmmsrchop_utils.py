"""
Utilities for the domtblop.py scripts.
"""

import os
import sys

from typing import (
    TextIO,
    Tuple,
)


def eprint(*args, **kwargs) -> None:
    """
    Print to stderr.

    Parameters:
        *args: The arguments to print.
        **kwargs: The keyword arguments to print.
    """
    print(*args, file=sys.stderr, **kwargs)


def read_input(file_path: str) -> TextIO:
    """
    Opens the input file or stdin if the file_path is "-".

    Parameters:
        file_path (str): The path to the input file.

    Returns:
        file_handle: A file handle to the opened file or stdin.

    Raises:
        FileNotFoundError: If the input file is not found.
    """

    if file_path == "-":
        return sys.stdin
    elif file_path.startswith("/dev/fd/"):
        fd = int(os.path.basename(file_path))
        return os.fdopen(fd, "r")
    else:
        if not os.path.isfile(file_path):
            eprint(f"Error: File not found: {file_path}")
            raise FileNotFoundError

        return open(file_path, "r")


def parse_seqkit2translate_fasta_header(header: str) -> Tuple[str, int, int, int]:
    """
    Parse the header of a seqkit2 translated FASTA file.
    The header format is expected to be:
        [parent_seq]_frame=[frame]_begin=[begin]_end=[end]

    Parameters:
        header (str): The header string to parse.

    Returns:
        Tuple[str, int, int, int]: A tuple containing:
            1. parent_seq (str): The parent sequence identifier.
            2. frame (int): The frame number.
            3. begin (int): The beginning position.
            4. end (int): The ending position.

    Raises:
        ValueError: If the header format is invalid.
    """
    if header.count("=") != 3:
        eprint(f"Error: Invalid header format: {header}")
        raise ValueError("Invalid header format")

    parent_seq_str, frame_str, begin_str, end_str = header.split("=")

    parent_seq = parent_seq_str[:-6]
    frame = int(frame_str[:-6])
    begin = int(begin_str[:-4])
    end = int(end_str)

    return parent_seq, frame, begin, end
