import os
from typing import Union
from Bio import SeqIO, SeqUtils
from abc import ABC, abstractmethod


class NucleotideNotFoundError(ValueError):
    """
    Error raised when there is incorrect nucleotide in corresponding sequence.
    """
    pass


class BiologicalSequence(ABC):
    @abstractmethod
    def __len__(self):
        """Return the length of the biological sequence."""
        pass

    def __getitem__(self, item):
        """Return the element or slice of the biological sequence."""
        pass

    def __str__(self):
        """Return the string representation of the biological sequence."""
        pass

    def is_correct(self):
        """Check if the alphabet of the sequence is correct."""
        pass


class NucleicAcidSequence(BiologicalSequence, ABC):
    def __init__(self, seq: str):
        self.seq = seq
        self.complement_pairs = {}

    def is_correct(self):
        """
        Check if seq dna/rna or not
        """
        for nucleotide in set(self.seq):
            if nucleotide.upper() not in self.complement_pairs.keys():
                raise NucleotideNotFoundError(f'Sequence contains invalid nucleotide: {nucleotide}')
        return None

    def complement(self):
        """
        Create complementary nucleic seq
        """
        self.is_correct()
        complement_seq = ''

        for nucleotide in self.seq:
            if nucleotide.isupper():
                complement_seq += self.complement_pairs[nucleotide]
            else:
                complement_seq += self.complement_pairs[nucleotide.upper()].lower()
        result = type(self)(complement_seq)
        return result

    def gc_content(self):
        self.is_correct()
        return SeqUtils.gc_fraction(self.seq)


class DNASequence(NucleicAcidSequence, BiologicalSequence, ABC):
    def __init__(self, seq):
        super().__init__(seq)
        self.complement_pairs = {'A': 'T',
                                 'T': 'A',
                                 'C': 'G',
                                 'G': 'C'}

    def transcribe(self):
        """
        Calculate rna transcript from dna seq

        :return:
        - str: rna transcript
        """
        if 'T' in self.seq or 't' in self.seq:
            return self.seq.replace('T', 'U').replace('t', 'u')
        return self.seq


class RNASequence(NucleicAcidSequence, BiologicalSequence, ABC):
    def __init__(self, seq):
        super().__init__(seq)
        self.complement_pairs = {'A': 'U',
                                 'U': 'A',
                                 'C': 'G',
                                 'G': 'C'}


def run_fastaq_filtering(input_path: str, output_filename: str = None, gc_bounds: Union[tuple, float, int] = (0, 100),
                         length_bounds: Union[tuple, float] = (0, 2 ** 32),
                         quality_threshold: int = 0) -> None:
    """
    Launch filtering fastaq seq using 3 adjustable cutoffs. Allowed intervals include cutoffs values.

    :param input_path: path to input fastaq file
    :param output_filename: name of output fastaq file
    :param gc_bounds: cutoff for GC content in percents. You can specify lower and upper limits (tuple with floats)
    or just upper limit (then pass float). Default = (0,100)
    :param length_bounds: cutoff fot length in nucleic bases. You can specify lower and upper limits (tuple with floats)
    or just upper limit (then pass float ot int). Default = (0,2**32)
    :param quality_threshold: cutoff for seq quality in phred33 scale. Default = 0.
    Reads with average score lower than this cutoff will be dropped.

    :return: None
    This function does not return anything. It saves the filtered FASTQ sequences
    in the specified output file in fastq_filtrator_results directory.
    """
    filtered_seqs = []

    for record in SeqIO.parse(input_path, 'fastq'):
        # create dict for collect filtration steps results
        filter_result = {}

        # filter based on GC content
        gc_percent = SeqUtils.gc_fraction(record) * 100
        if isinstance(gc_bounds, tuple):
            filter_result['gc_content'] = (gc_bounds[0] <= gc_percent <= gc_bounds[1])
        else:
            filter_result['gc_content'] = (gc_percent <= gc_bounds)

        # filter based on length
        length = len(record.seq)
        if isinstance(length_bounds, tuple):
            filter_result['length'] = (length_bounds[0] <= length <= length_bounds[1])
        else:
            filter_result['length'] = (length <= length_bounds)

        # filter based on quality
        quality = sum(record.letter_annotations['phred_quality']) / len(record)
        filter_result['quality'] = (quality >= quality_threshold)

        # total filtering: take structures that passed all cutoffs and write them to the output file
        if filter_result['gc_content'] and filter_result['length'] and filter_result['quality']:
            filtered_seqs.append(record)

    os.makedirs('fastq_filtrator_results', exist_ok=True)

    SeqIO.write(filtered_seqs, os.path.join('fastq_filtrator_results', output_filename), 'fastq')


if __name__ == "__main__":
    path_to_fasta = '/home/polina/bioinf/python_programming/HW6_Files/example_data/test.fastq'
    run_fastaq_filtering(path_to_fasta, 'test.fastq', gc_bounds=(3, 50), quality_threshold=20)
