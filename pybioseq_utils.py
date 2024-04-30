import datetime
import os
from abc import ABC, abstractmethod
from typing import Union

from Bio import SeqIO, SeqUtils


class NucleotideNotFoundError(ValueError):
    """
    Error raised when there is incorrect nucleotide in corresponding sequence.
    """

    pass


class AminoAcidNotFoundError(ValueError):
    """
    Error raised when there is incorrect aminoacid in corresponding sequence.
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

    def check_seq_correctness(self):
        """Describe protein seq biological sequence correctness."""
        pass


class NucleicAcidSequence(BiologicalSequence, ABC):
    """
    Represents a nucleic acid sequence.

    Inherits from:
        BiologicalSequence: An abstract base class defining the interface for biological sequences.

    Attributes:
        seq (str): The nucleic acid sequence.
        complement_pairs (dict): A dictionary mapping nucleotides to their complements.

    Methods:
        __len__(): Returns the length of the nucleic acid sequence.
        __getitem__(item): Returns the element or slice of the nucleic acid sequence.
        __str__(): Returns the string representation of the nucleic acid sequence.
        check_seq_correctness(): Check if the nucleic acid sequence contains valid nucleotides and raise error if not.
        complement(): Creates a complementary nucleic acid sequence.
        gc_content(): Calculates the GC content of the sequence.
    """

    complement_pairs = {}

    def __init__(self, seq: str):
        self.seq = seq

    def __len__(self):
        return len(self.seq)

    def __getitem__(self, index):
        return self.seq[index]

    def __str__(self):
        return str(self.seq)

    def check_seq_correctness(self):
        """
        Describe nucleic seq biological sequence correctness and raise error in invalid letters were provided.
        """
        for nucleotide in set(self.seq):
            if nucleotide.upper() not in self.complement_pairs.keys():
                raise NucleotideNotFoundError(
                    f"Sequence contains invalid nucleotide: {nucleotide}"
                )
        return "Your sequence is correct."

    def complement(self):
        """
        Create complementary nucleic acid sequence.
        """
        self.check_seq_correctness()
        complement_seq = ""

        for nucleotide in self.seq:
            if nucleotide.isupper():
                complement_seq += self.complement_pairs[nucleotide]
            else:
                complement_seq += self.complement_pairs[nucleotide.upper()].lower()
        resulting_seq = type(self)(complement_seq)
        return resulting_seq

    def gc_content(self):
        """
        Calculate the GC content of the sequence as fraction (from 0 to 1).
        """
        self.check_seq_correctness()
        return SeqUtils.GC(self.seq)


class RNASequence(NucleicAcidSequence, BiologicalSequence, ABC):
    """
    Represents an RNA sequence.

    Inherits from:
        NucleicAcidSequence: A generic nucleic acid sequence.

    Methods:
        __init__(): Initializes an RNA sequence.
    """

    complement_pairs = {"A": "U", "U": "A", "C": "G", "G": "C"}

    def __init__(self, seq):
        super().__init__(seq)


class DNASequence(NucleicAcidSequence, BiologicalSequence, ABC):
    """
    Represents a DNA sequence.

    Inherits from:
        NucleicAcidSequence: A generic nucleic acid sequence.

    Methods:
        __init__(): Initializes a DNA sequence.
        transcribe(): Transcribes the DNA sequence into RNA.
    """

    complement_pairs = {"A": "T", "T": "A", "C": "G", "G": "C"}

    def __init__(self, seq):
        super().__init__(seq)

    def transcribe(self):
        """
        Calculate rna transcript from dna seq
        """
        if "T" in self.seq or "t" in self.seq:
            transcript = type(self)(self.seq.replace("T", "U").replace("t", "u"))
        else:
            transcript = type(self)(self.seq)

        return transcript


class AminoAcidSequence(BiologicalSequence, ABC):
    """
    Represents an amino acid sequence.

    Inherits from:
        BiologicalSequence: An abstract base class defining the interface for biological sequences.

    Attributes:
        seq (str): The amino acid sequence.
        code (int): The encoding format for the sequence (1 for one-letter and 3 for three-letter).

    Methods:
    __len__(): Returns the length of the amino acid sequence. __getitem__(index): Returns the amino acid at
    the specified index.
    __str__(): Returns the string representation of the amino acid sequence.
    check_seq_correctness(): Describe the amino acid sequence correctness and raise error in it contains invalid
    letters.
    count_aa(): Counts the occurrences of each amino acid in the sequence.
    """

    def __init__(self, seq: str, code: int = 1):
        self.code = code
        self.residues_encoding = {
            "ALA": "A",
            "ARG": "R",
            "ASN": "N",
            "ASP": "D",
            "CYS": "C",
            "GLN": "Q",
            "GLU": "E",
            "GLY": "G",
            "HIS": "H",
            "ILE": "I",
            "LEU": "L",
            "LYS": "K",
            "MET": "M",
            "PHE": "F",
            "PRO": "P",
            "SER": "S",
            "THR": "T",
            "TRP": "W",
            "TYR": "Y",
            "VAL": "V",
        }

        if self.code == 1:
            self.seq = [seq[val : val + 1] for val in range(0, len(seq), 1)]
            self.residues_names = self.residues_encoding.values()
        elif self.code == 3:
            self.seq = [seq[val : val + 3] for val in range(0, len(seq), 3)]
            self.residues_names = self.residues_encoding.keys()
        else:
            raise ValueError(
                f"Invalid specified encoding: {self.code}. Available values: 1, 3."
            )

    def __len__(self):
        return len(self.seq)

    def __getitem__(self, index):
        return self.seq[index]

    def __str__(self):
        return "".join(self.seq)

    def check_seq_correctness(self):
        """
        Describe protein seq correctness
        :return: seq description (str)
        """
        for residue in set(self.seq):
            if residue.upper() not in self.residues_names:
                raise AminoAcidNotFoundError(
                    f"Sequence contains invalid aminoacid: {residue}"
                )
        return "Your sequence is correct."

    def count_aa(self):
        """
        Count entry of each residue type in your seq. Get description of amino acid composition in dict format.
        :return: each residue type and its amount in current seq (dict)
        """
        residue_count = {}
        for residue in set(self.seq):
            residue_entry = self.seq.count(residue)
            residue_count[residue] = residue_entry
        return residue_count


def filter_fastq(
    input_path: str,
    output_filename: str = None,
    gc_bounds: Union[tuple, float, int] = (0, 100),
    length_bounds: Union[tuple, float] = (0, 2**32),
    quality_threshold: int = 0,
) -> None:
    """
    Launch filtering fastq seq using 3 adjustable cutoffs. Allowed intervals include cutoffs values.

    :param input_path: path to input fastq file
    :param output_filename: name of output fastq file
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

    for record in SeqIO.parse(input_path, "fastq"):
        # create dict for collect filtration steps results
        filter_result = {}

        # filter based on GC content
        gc_percent = SeqUtils.gc_fraction(record) * 100
        if isinstance(gc_bounds, tuple):
            filter_result["gc_content"] = gc_bounds[0] <= gc_percent <= gc_bounds[1]
        else:
            filter_result["gc_content"] = gc_percent <= gc_bounds

        # filter based on length
        length = len(record.seq)
        if isinstance(length_bounds, tuple):
            filter_result["length"] = length_bounds[0] <= length <= length_bounds[1]
        else:
            filter_result["length"] = length <= length_bounds

        # filter based on quality
        quality = sum(record.letter_annotations["phred_quality"]) / len(record)
        filter_result["quality"] = quality >= quality_threshold

        # total filtering: take structures that passed all cutoffs and write them to the output file
        if (
            filter_result["gc_content"]
            and filter_result["length"]
            and filter_result["quality"]
        ):
            filtered_seqs.append(record)

    os.makedirs("fastq_filtrator_results", exist_ok=True)

    # create filename if it is not given
    if output_filename is None:
        f"fasta_filtrator_{datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}_out"

    SeqIO.write(
        filtered_seqs, os.path.join("fastq_filtrator_results", output_filename), "fastq"
    )
