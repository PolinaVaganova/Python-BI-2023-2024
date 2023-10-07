from typing import Union, Sequence, List

DNA_DNA_PAIRS = {'A': 'T',
                 'T': 'A',
                 'C': 'G',
                 'G': 'C'}

RNA_DNA_PAIRS = {'A': 'T',
                 'U': 'A',
                 'C': 'G',
                 'G': 'C'}

RNA_RNA_PAIRS = {'A': 'U',
                 'U': 'A',
                 'C': 'G',
                 'G': 'C'}

NUCLEOTIDE_NAMES = ['A', 'T', 'C', 'G', 'U']


# function for getting a rna transcript, based on dna sequence
def transcribe(seq: str) -> str:
    """
    Calculate rna transcript seq from dna seq

    :param seq:
    - seq (str): dna seq

    :return:
    - str: rna transcript seq
    """
    if 'T' in seq or 't' in seq:
        transcript = seq.replace('T', 'U').replace('t', 'u')
    elif 'U' in seq.upper():
        return 'Passed rna to function, dna excepted! Skip it!'
    else:
        transcript = seq
    return transcript


# function for getting a dna, based on rna transcript
def reverse_transcribe(seq: str) -> str:
    """
    Calculate complementary dna seq from rna seq

    :param seq:
    - seq (str): rna seq

    :return:
    - str: complementary dna seq
    """
    transcript = ''
    for nucleotide in seq:
        if nucleotide.upper() == 'T':
            return 'Passed dna to function, rna excepted! Skip it!'
        elif nucleotide.isupper():
            transcript += RNA_DNA_PAIRS[nucleotide]
        else:
            transcript += RNA_DNA_PAIRS[nucleotide.upper()].lower()
    return transcript


def reverse(seq: str) -> str:
    """
    Get reverse seq for dna or rna

    :param seq:
    - seq (str): dna or rna seq

    :return:
    - str: reverse seq
    """
    pass


# function for getting a complementary sequence for dna or rna
def complement(seq: str) -> str:
    """
    Create complementary seq for dna or rna (if unclear will be treated as dna)

    :param seq:
    - seq (str): dna or rna seq

    :return:
    - str: complementary seq
    """
    pass


# function for getting a complementary sequence for dna or rna in reverse format
def reverse_complement(seq: str) -> str:
    """
    Create complementary seq for dna or rna in reverse format (if unclear will be treated as dna)

    :param seq:
    - seq (str): dna or rna seq

    :return:
    - str: complementary seq in reverse format
    """
    pass


# function for count nucleotides in seq
def count_nucleotides(seq: str) -> str:
    """
    Characterize percentages of each type of nucleotide

    :param seq:
    - seq (str): dna or rna seq

    :return:
    - str: each nucleotide with it percentage weight
    """
    pass


def make_triplets(seq: str) -> list:
    """
    Split your seq into triplets

    :param:
    - seq (str): dna or rna seq

    :return:
    - List[str]: created triplets
    """
    pass


# function for check a sequence correctness
def is_dna_or_rna(seq: str) -> bool:
    """
    Check if seq dna/rna or not

    :param:
    - seq (str): dna or rna seq

    :return:
    - bool: the result of the check
    """
    pass
