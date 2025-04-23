#!/usr/bin/env python

"""
Build a test dataset for the 'tgscan.nf' and 'exschuf.nf' pipelines.
"""

from dataclasses import dataclass
import sys
from typing import (
    List, Tuple, Dict
)


OUT_FASTA = "test_dataset.fasta"
OUT_GFF = "test_dataset.gff"


@dataclass
class SyntheticDnaSequence:
    """
    This class defines a synthetic test dataset for LRR analysis.
    """

    # Dataset name
    source: str = "synthetic"
    sequence: str = "TTTCCTCTCGAGTCCCGTCCAGTTGAGCGTATCACTCCCAGTGTACTAGCAAGCCGAGAAGGCTGTGCTTGGAGTCAATCGGATGTAGGATGGTCTCCAGACACCGGGCCACCACTCTTCACGCCTAAAGCATAAACGTCGAGCAGTCATGAAAGTCTTAGTACCGGACGTGCCGTTTCACTGCGAATATTACCTGAAGCTGTACCGTTATTGCGGAGCAAAGATGCAGTGCTGCTCTTATCATATTTGTATTGACGACAGCCGCCTTCGCGGTTTCCTCAGACACTTAAGAATAAGGGCTTATTGTAGGCAGAGGCACGCCCTTTTAGTGGCTGCGGCAAAATATCTTCGGATCCCCTTGTCTAACCAAATTAATCGAATTCTCTCATTTAAGACCCTAATATGTCATCATTAGTGTTTAAATGCCACCCCGAAAATACCGCCTAGAAATGTCTATGATTGGTCCACTAAAGTTGATTAAAACGACTGCTAAATCCGCGTGATAGGGCATTTGAAGTTTAATTTTGTATCGCAAGGTACTCCCGATCTTAATGGATGGCCGGAAGTGGTACGGATGCAATAAGCGCGGGTGAGAGGGTAATTAGGCGCGTTCACCTACGCTACGCTAACGGGCGATTCTATAAGAATGCACATTGCGTCGATTCATAAGATGTCTCGACCGCATGCGCAACTTGTGAAGTGTCTACTATCCCTAAGCGCATATCTCGCACAGTAACCCCCGAATATGTCGGCATCTGATGTTACCCGGGTTGAGTTAGTGTTGAGCTCACGGAACTTATTGTATGAGTAGAGATTTGTAAGAGCTGTTAGTTAGCTCGCTCAGCTAATAGTTGCCCACACAACGTCAAAATTAGAGAACGGTCGTAACATTATCGGTGGTTCTCTAACTACTATCAGTACCCACGACTCGACTCTGCCGCAGCTACGTATCGCCTGAAAGCCAGTCAGCGTTAAGGAGTGCTCTGACCAGGACAACACG"

@dataclass
class LrrMorif:
    """
    This class defines a synthetic test dataset for LRR analysis.
    """

    # Dataset name
    source: str = "Ec-03_001790.1"
    lrr_morif_protein: str = "PILKELGALTKLTSLFLRSNKLTG"
    lrr_morif_dna: str = "CCGATTCTGAAAGAACTGGGCGCGCTGACCAAACTGACCAGCCTGTTTCTGCGCAGCAACAAACTGACCGGC"
    bp_length: int = 72

@dataclass
class NonLrrMorif:
    """
    This class defines a synthetic test dataset for LRR analysis.
    """

    # Dataset name
    source: str = "LexA_Binding_Domain"
    non_lrr_morif_protein: str = "MKALTARQQEVFDLIRDHISQTGMPPTRAEIAQRLGFRSPNAAEEHLKALARKGVIEIVS"
    non_lrr_morif_dna: str = "ATGAAAGCGCTGACCGCGCGCCAGCAGGAAGTGTTTGATCTGATTCGCGATCATATTAGCCAGACCGGCATGCCGCCGACCCGCGCGGAAATTGCGCAGCGCCTGGGCTTTCGCAGCCCGAACGCGGCGGAAGAACATCTGAAAGCGCTGGCGCGCAAAGGCGTGATTGAAATTGTGAGC"
    bp_length: int = 180

@dataclass
class Fasta:
    header: str
    sequence: str
    _wrap: int = 80

    def wrap_sequence(self):
        """
        Wrap the sequence to the specified line length.
        """
        wrapped_sequence = []
        for i in range(0, len(self.sequence), self._wrap):
            wrapped_sequence.append(self.sequence[i:i + self._wrap])
        return "\n".join(wrapped_sequence)

    def __str__(self):
        return f">{self.header}\n{self.wrap_sequence()}"

@dataclass
class GffFeature:
    seqid: str
    source: str
    type_: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attributes: Dict[str, str]

    def __str__(self):
        attributes_str = ";".join(f"{k}={v}" for k, v in self.attributes.items())
        return f"{self.seqid}\t{self.source}\t{self.type_}\t{self.start}\t{self.end}\t{self.score}\t{self.strand}\t{self.phase}\t{attributes_str}"


def eprint(*args, **kwargs):
    """
    Print to stderr
    """
    print(*args, file=sys.stderr, **kwargs)

def mk_chr_1() -> Tuple[Fasta, List[GffFeature]]:
    """
    Chromosome 1 contains two genes:
        - Gene 1: Simple LRR gene, contains 1 single LRR exon
        - Gene 2: Simple non-LRR gene, contains 1 single non-LRR exon

    Total chromosome length: 3000 bp
        - Filler 1: 500 bp
        - Gene 1: 1000 bp
        - Gene 2: 1000 bp
        - Filler 2: 500 bp
    """
    gff_features = []

    synthetic_dna = SyntheticDnaSequence()
    lrr_morif = LrrMorif()
    non_lrr_morif = NonLrrMorif()

    # Filler 1
    filler_1 = synthetic_dna.sequence[:500]

    # Gene 1
    gene_1_seq = [
            synthetic_dna.sequence[:499],
            lrr_morif.lrr_morif_dna,
            synthetic_dna.sequence[:428],
            ]
    gene_1_seq = "".join(gene_1_seq)

    gene1_gff_mrna = GffFeature(
        seqid="chr1",
        source="synthetic",
        type_="mRNA",
        start=501,
        end=1500,
        score=".",
        strand="+",
        phase="0",
        attributes={"ID": "gene1"}
    )
    gff_features.append(gene1_gff_mrna)
    gene1_gff_cds = GffFeature(
        seqid="chr1",
        source="synthetic",
        type_="CDS",
        start=1000,
        end=1072,
        score=".",
        strand="+",
        phase="0",
        attributes={"Parent": "gene1"}
    )
    gff_features.append(gene1_gff_cds)


    # Gene 2
    gene_2_seq = [
            synthetic_dna.sequence[:500],
            non_lrr_morif.non_lrr_morif_dna,
            synthetic_dna.sequence[:320],
            ]
    gene_2_seq = "".join(gene_2_seq)
    gene2_gff_mrna = GffFeature(
        seqid="chr1",
        source="synthetic",
        type_="mRNA",
        start=1501,
        end=2500,
        score=".",
        strand="+",
        phase="0",
        attributes={"ID": "gene2"}
    )
    gff_features.append(gene2_gff_mrna)
    gene2_gff_cds = GffFeature(
        seqid="chr1",
        source="synthetic",
        type_="CDS",
        start=2000,
        end=2180,
        score=".",
        strand="+",
        phase="1",
        attributes={"Parent": "gene2"}
    )

    # Filler 2
    filler_2 = synthetic_dna.sequence[:500]

    # Create the chromosome sequence
    chromosome_seq = "".join([filler_1, gene_1_seq, filler_2, gene_2_seq])
    fasta = Fasta(
        header="chr1",
        sequence=chromosome_seq
    )

    return fasta, gff_features


def mk_test_dataset() -> Tuple[Fasta, List[GffFeature]]:
    """
    Create a test dataset for the 'tgscan.nf' and 'exschuf.nf' pipelines.
    """
    fasta, gff_features = mk_chr_1()

    # Write the FASTA file
    with open(OUT_FASTA, "w") as fasta_file:
        fasta_file.write(str(fasta))

    # Write the GFF file
    with open(OUT_GFF, "w") as gff_file:
        for feature in gff_features:
            gff_file.write(str(feature) + "\n")

    return fasta, gff_features


if __name__ == "__main__":
    # Create the test dataset
    fasta, gff_features = mk_test_dataset()

    # Print the results
    eprint(f"Printed fasta file: {OUT_FASTA}")
    eprint(f"Printed gff file: {OUT_GFF}")
