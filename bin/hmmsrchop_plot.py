#!/usr/bin/env python

"""
Generate informative plots of the tgscan.nf pipeline output.

Usage:
    hmmsrchop.py plot <serialized HMMER3 output> [options]

Required Arguments:
    <serialized HMMER3 output>  The serialized HMMER3 output file to convert to GFF format.

Optional Arguments:
    --seq-eval-threshold <float>        Per sequence E-value threshold to use for filtering the results. Default: null
    --seq-bitscore-threshold <float>    Per sequence bitscore threshold to use for filtering the results. Default: null
    --align-eval-threshold <float>      Per profile alignment E-value threshold to use for filtering the results. Default: null
    --align-bitscore-threshold <float>  Per profile allignment bitscore threshold to use for filtering the results. Default: null
"""

import argparse
from collections import defaultdict
from itertools import chain
import os
import sys
from typing import List, Dict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from hmmsrchop_parse import TargetSequence
from hmmsrchop_utils import (
    eprint,
    read_input,
)

MIN_EVAL = 1e-30
MAX_BITSCORE = 100.0


# Set the default style for the plots
sns.set_theme()

# Set up the figure sizes
pt = 1.0 / 72.27  # 72.27 points per inch
two_columns_figure = 539.0 * pt
golden_ratio = (1 + 5**0.5) / 2

my_width = two_columns_figure

HORIZONTAL_FIGSIZE = (my_width, my_width / golden_ratio)
VERTICAL_FIGSIZE = (my_width * golden_ratio, my_width)
XXLARGE_HORIZ_FIGSIZE = (my_width * 5, my_width * 5 / golden_ratio)
XXLARGE_VERT_FIGSIZE = (my_width * 5 * golden_ratio, my_width * 5)


def evalue_distr(
    evalues: List[float],
    filepath: str,
    title: str,
    threshold: float | None = None,
) -> None:
    # Check type of evalues is float
    if not all(isinstance(i, float) for i in evalues):
        eprint("Error: E-values must be a list of floats")
        raise ValueError("evalues must be a list of floats")

    df = pd.DataFrame({"evalue": evalues})
    df["evalue"] = df["evalue"].replace(0, MIN_EVAL)
    df["-log10_evalue"] = -np.log10(df["evalue"])

    f, ax = plt.subplots(figsize=HORIZONTAL_FIGSIZE)

    sns.set_color_codes("muted")

    sns.histplot(
        x="-log10_evalue",
        data=df,
        ax=ax,
        kde=True,
    )

    ax.set_title(title)
    ax.set_ylabel("Count")
    ax.set_xlabel("-log10(E-value)")

    # Set the x-axis limits to be between 0 and 100
    ax.set_xlim(0.0, -np.log10(MIN_EVAL))

    if threshold:
        ax.axvline(
            x=float(-np.log10(threshold)),
            color="red",
            linestyle="--",
            label=f"Threshold: <={threshold}",
        )
        ax.legend()

    sns.despine()

    plt.savefig(f"{filepath}.jpeg", dpi=300, bbox_inches="tight", format="jpeg")
    plt.clf()
    plt.close()


def evalue_distr_per_profile(
    evalues: Dict[str, List[float]],
    filepath: str,
    threshold: float | None = None,
) -> None:
    # helper function to annotate the plot with the threshold
    def annotate_plot(*args, **kwargs):
        for ax in g.axes.flat:
            x_label_loc = -np.log10(threshold)
            # Draw a vertical line
            ax.axvline(
                x=x_label_loc,
                color="r",
                linestyle="--",
                label=f"Threshold: <={threshold}",
            )
            # Add a label
            ax.text(
                x=x_label_loc,
                y=ax.get_ylim()[1] * 0.9,
                s=f"Threshold: <={threshold}",
                color="r",
                fontsize=10,
                horizontalalignment="left",
                verticalalignment="center",
            )

    data = {"profile": [], "evalue": []}
    for profile, value in evalues.items():
        data["profile"].extend([profile] * len(value))
        data["evalue"].extend(value)

    df = pd.DataFrame(data)
    df["evalue"] = df["evalue"].replace(0, MIN_EVAL)

    df["-log10_evalue"] = -np.log10(df["evalue"])

    unique_profiles = df["profile"].unique().tolist()
    if len(unique_profiles) > 10:
        col_wrap = 3
    elif len(unique_profiles) < 17:
        col_wrap = 4
    elif len(unique_profiles) < 26:
        col_wrap = 5
    elif len(unique_profiles) < 37:
        col_wrap = 6
    else:
        raise ValueError("Too many profiles to plot")

    f, ax = plt.subplots(figsize=HORIZONTAL_FIGSIZE)

    sns.set_color_codes("muted")
    g = sns.displot(
        x="-log10_evalue",
        data=df,
        kde=True,
        col="profile",
        col_wrap=col_wrap,
        facet_kws={"sharey": False, "sharex": False},
    )
    g.set(xlim=(0.0, -np.log10(MIN_EVAL)))

    if threshold:
        g.map(annotate_plot)

    sns.despine()
    plt.savefig(f"{filepath}.jpeg", dpi=300, bbox_inches="tight", format="jpeg")
    plt.clf()
    plt.close()


def bitscore_distr(
    bitscores: List[float],
    filepath: str,
    title: str,
    threshold: float | None = None,
) -> None:
    # Check type of bitscores is float
    if not all(isinstance(i, float) for i in bitscores):
        eprint("Error: Bitscores must be a list of floats")
        raise ValueError("bitscores must be a list of floats")

    df = pd.DataFrame({"bitscore": bitscores})

    f, ax = plt.subplots(figsize=HORIZONTAL_FIGSIZE)

    sns.set_color_codes("muted")

    sns.histplot(
        x="bitscore",
        data=df,
        ax=ax,
        kde=True,
    )

    ax.set_title(title)
    ax.set_ylabel("Count")
    ax.set_xlabel("Bitscore")

    if threshold:
        ax.axvline(
            x=threshold,
            color="red",
            linestyle="--",
            label=f"Threshold: >={threshold}",
        )
        ax.legend()

    sns.despine()

    plt.savefig(f"{filepath}.jpeg", dpi=300, bbox_inches="tight", format="jpeg")
    plt.clf()
    plt.close()


def setup_argparse() -> argparse.ArgumentParser:
    """
    Sets up the argparse instance for command-line arguments.

    Returns:
    argparse.ArgumentParser: Configured ArgumentParser instance.
    """
    fmt = lambda prog: argparse.RawTextHelpFormatter(prog)

    parser = argparse.ArgumentParser(
        add_help=False,
        formatter_class=fmt,
        description=__doc__,
    )

    parser.add_argument(
        "input",
        type=str,
        help="The serialized HMMER3 output file to convert to GFF format.",
        default="-",
    )

    parser.add_argument(
        "--seq-eval-threshold",
        type=float,
        help="Per sequence E-value threshold to use for filtering the results. Default: null",
        default=None,
    )

    parser.add_argument(
        "--seq-bitscore-threshold",
        type=float,
        help="Per sequence bitscore threshold to use for filtering the results. Default: null",
        default=None,
    )

    parser.add_argument(
        "--align-eval-threshold",
        type=float,
        help="Per profile alignment E-value threshold to use for filtering the results. Default: null",
        default=None,
    )

    parser.add_argument(
        "--align-bitscore-threshold",
        type=float,
        help="Per profile allignment bitscore threshold to use for filtering the results. Default: null",
        default=None,
    )

    # General Options
    parser.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS)

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

    perseq_evalues = defaultdict(list)
    perseq_bitscores = defaultdict(list)
    peraln_evalues = defaultdict(list)
    peraln_bitscores = defaultdict(list)
    for target_sequence in target_sequences:
        for profile in target_sequence.hits:
            perseq_evalues[profile.name].append(profile.evalue)
            perseq_bitscores[profile.name].append(profile.bitscore)

            for alignment in profile.alignments:
                peraln_evalues[profile.name].append(alignment.conditional_evalue)
                peraln_bitscores[profile.name].append(alignment.bitscore)

    if len(list(chain.from_iterable(perseq_evalues.values()))) < 10:
        eprint("Error: Not enough data to plot.")
        sys.exit(0)

    if len(list(chain.from_iterable(peraln_evalues.values()))) < 10:
        eprint("Error: Not enough data to plot.")
        sys.exit(0)

    # Plot the distribution of E-values for all profiles
    evalue_distr(
        list(chain.from_iterable(perseq_evalues.values())),
        f"all_profiles_per_sequence_evalue_distr",
        "Distribution of E-values for all profiles",
        config.seq_eval_threshold,
    )

    # Plot the distribution of bitscores for all profiles
    bitscore_distr(
        list(chain.from_iterable(perseq_bitscores.values())),
        f"all_profiles_per_sequence_bitscore_distr",
        "Distribution of bitscores for all profiles",
        config.seq_bitscore_threshold,
    )

    # Plot the distribution of E-values for all alignments
    evalue_distr(
        list(chain.from_iterable(peraln_evalues.values())),
        f"all_alignments_per_sequence_evalue_distr",
        "Distribution of E-values for all alignments",
        config.align_eval_threshold,
    )

    # Plot the distribution of bitscores for all alignments
    bitscore_distr(
        list(chain.from_iterable(peraln_bitscores.values())),
        f"all_alignments_per_sequence_bitscore_distr",
        "Distribution of bitscores for all alignments",
        config.align_bitscore_threshold,
    )

    evalue_distr_per_profile(
        perseq_evalues,
        f"per_profile_evalue_distr",
        config.seq_eval_threshold,
        )
