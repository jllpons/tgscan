"""
Utilities for the domtblop.py scripts.
"""

import os
import sys

from typing import (
    TextIO,
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
