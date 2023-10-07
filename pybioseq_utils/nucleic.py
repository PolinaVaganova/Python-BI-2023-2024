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
    return seq[::-1]


# function for getting a complementary sequence for dna or rna
def complement(seq: str) -> str:
    """
    Create complementary seq for dna or rna (if unclear will be treated as dna)

    :param seq:
    - seq (str): dna or rna seq

    :return:
    - str: complementary seq
    """
    complement_seq = ''
    if 'U' in seq:
        for nucleotide in seq:
            if nucleotide.isupper():
                complement_seq += RNA_RNA_PAIRS[nucleotide]
            else:
                complement_seq += RNA_RNA_PAIRS[nucleotide.upper()].lower()
    else:
        for nucleotide in seq:
            if nucleotide.isupper():
                complement_seq += DNA_DNA_PAIRS[nucleotide]
            else:
                complement_seq += DNA_DNA_PAIRS[nucleotide.upper()].lower()
    return complement_seq


# function for getting a complementary sequence for dna or rna in reverse format
def reverse_complement(seq: str) -> str:
    """
    Create complementary seq for dna or rna in reverse format (if unclear will be treated as dna)

    :param seq:
    - seq (str): dna or rna seq

    :return:
    - str: complementary seq in reverse format
    """
    return reverse(complement(seq))


# function for count nucleotides in seq
def count_nucleotides(seq: str) -> dict:
    """
    Find amount of each type of nucleotide in seq

    :param seq:
    - seq (str): dna or rna seq

    :return:
    - dict: each nucleotide with it entry number
    """
    content = {'A': 0,
               'T': 0,
               'C': 0,
               'G': 0,
               'U': 0}
    for nucleotide in seq:
        content[nucleotide.upper()] += 1
    return content


def make_triplets(seq: str) -> List[str]:
    """
    Split your seq into triplets

    :param:
    - seq (str): dna or rna seq

    :return:
    - List[str]: created triplets
    """
    triplets = []
    for ind in range(0, len(seq) - 2):
        triplet = seq[ind:ind + 4]
        triplets.append(triplet)
    return triplets


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
