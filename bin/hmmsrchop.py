#!/usr/bin/env python

"""
Parse and serialize hmmsearch output file and perform different operations

Usage: hmmsrchop.py <operation> [options]

    Parsing and Serialization:
        parse       Parse and serialize hmmsearch into target sequence-centric JSON structures.

    Operations on Serialized Domtblout Target Results:
        togff       Convert serialized domain hits to GFF format.
        tobed       Convert serialized domain hits to BED format.
        topkl       Convert serialized domain hits to Pickle format.

    Utilities:
        plot        Generate some plots in `$PWD/plots/` directory.
        help        Show this help message and exit.

Options:
    -h, --help      Show this help message and exit.
    -v, --version   Show version and exit.
"""


import sys


__version__ = "0.1.0"


def main():
    argn = len(sys.argv)
    if argn == 1:
        print(__doc__)
        sys.exit(1)

    cmd = sys.argv[1]


    for i in ["-v", "--version"]:
        if i in sys.argv:
            print(f"domtblop.py v{__version__}")
            sys.exit(1)

    if cmd == "-h" or cmd == "--help" or cmd == "help":
        print(__doc__)
        sys.exit(1)

    elif cmd == "parse":
        from hmmsrchop_parse import run

        run(sys.argv[2:])

    elif cmd == "togff":
        from domtblop_togff import run

        run(sys.argv[2:])

    elif cmd == "tobed":
        from domtblop_tobed import run

        run(sys.argv[2:])

    elif cmd == "plot":
        from domtblop_plot import run

        run(sys.argv[2:])

    else:
        print(__doc__)
        print(f"[domtblop.py] Error: Unknown operation '{cmd}'", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
