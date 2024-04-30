import datetime
import io
import os
import re
import sys
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List, Optional
from typing import Union

import requests
from Bio import SeqIO, SeqUtils
from bs4 import BeautifulSoup


# Utils for biological sequence analysis


# Custom errors for sequences
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


# Abstract cass for sequences
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


# Class for nucleic acids sequences
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


# Class for amino acids sequences
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


# Function for fastaq files filtering
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


# Telegram bot for logging
def telegram_logger(chat_id: int):
    """
    Decorator function for logging messages to a Telegram chat.

    Args:
        chat_id (int): The ID of the Telegram chat to which the messages will be sent.

    Returns:
        function: Decorator function.
    """

    def decorator(func):
        def inner_func(*args, **kwargs):
            # Telegram bot token
            bot_token = os.getenv("TELEGRAM_BOT_TOKEN")
            if bot_token is None:
                print(
                    "Telegram bot token is not set. Please set the TELEGRAM_BOT_TOKEN environment variable."
                )
                return

            # Create a buffer to capture stdout and stderr
            buffer_stdout = io.StringIO()
            buffer_stderr = io.StringIO()
            sys.stdout = buffer_stdout
            sys.stderr = buffer_stderr

            start_time = datetime.datetime.now()

            func_name = func.__name__

            try:
                func(*args, **kwargs)
            except Exception as error:
                end_time = datetime.datetime.now()
                # Construct error message
                message = (
                    f"ERROR in <code>{func_name}</code>:\n"
                    f"<code>{repr(error)}</code>\n"
                    f"Process execution time: {end_time - start_time}"
                )
                # Send message to Telegram
                _send_telegram_message(
                    bot_token, chat_id, message, buffer_stdout, buffer_stderr, func_name
                )
                # Propagate the exception
                raise

            else:
                end_time = datetime.datetime.now()
                # Construct success message
                message = (
                    f"Process <code>{func.__name__}</code> completed\n"
                    f"Process execution time: {end_time - start_time}"
                )
                # Send message to Telegram
                _send_telegram_message(
                    bot_token, chat_id, message, buffer_stdout, buffer_stderr, func_name
                )

            finally:
                # Restore stdout and stderr
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__

        return inner_func

    return decorator


def _send_telegram_message(
    bot_token, chat_id, message, buffer_stdout, buffer_stderr, func_name
):
    """
    Sends a message and log document to a Telegram chat.

    Args:
        bot_token (str): The Telegram bot token.
        chat_id (int): The ID of the Telegram chat to which the message and log document will be sent.
        message (str): The message content to be sent.
        buffer_stdout (io.StringIO): The buffer containing stdout log information.
        buffer_stderr (io.StringIO): The buffer containing stderr log information.
        func_name (str): The name of the function associated with the message and log.

    Returns:
        None
    """
    # Send message to Telegram
    requests.post(
        f"https://api.telegram.org/bot{bot_token}/sendMessage",
        data={"chat_id": chat_id, "text": message, "parse_mode": "HTML"},
    )

    # Send stdout log document to Telegram
    stdout_log_name = f"log_{func_name}_stdout_{datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.log"
    stdout_log = (stdout_log_name, buffer_stdout.getvalue())
    requests.post(
        f"https://api.telegram.org/bot{bot_token}/sendDocument",
        data={"chat_id": chat_id},
        files={"document": stdout_log},
    )

    # Send stderr log document to Telegram
    stderr_log_name = f"log_{func_name}_stderr_{datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.log"
    stderr_log = (stderr_log_name, buffer_stderr.getvalue())
    requests.post(
        f"https://api.telegram.org/bot{bot_token}/sendDocument",
        data={"chat_id": chat_id},
        files={"document": stderr_log},
    )


# Bioinformatics API for GenScan requests
@dataclass
class GenscanOutput:
    """
    Dataclass representing the output of the Genscan analysis.

    Attributes:
        status (str): The status code of the HTTP response.
        cds_list (List[str]): List of predicted protein sequences.
        intron_list (List[List[Any]]): List of intron coordinates.
        exon_list (List[List[Any]]): List of exon coordinates.
    """

    status: str
    cds_list: List[str]
    intron_list: List[List[str]]
    exon_list: List[List[str]]


def extract_cds_list(results_text: str) -> List[str]:
    """
    Extracts predicted protein sequences from the Genscan results text.

    Args:
        results_text (str): Text containing Genscan results.

    Returns:
        List[str]: List of predicted protein sequences.
    """
    cds_list = re.findall(r"^>.+\n([\w\s\n]+)", results_text, flags=re.MULTILINE)
    cds_list = [sequence.replace("\n", "") for sequence in cds_list]

    return cds_list


def extract_exon_list(results_text) -> List[List[str]]:
    """
    Extracts exon coordinates from the Genscan results text.

    Args:
        results_text (str): Text containing Genscan results.

    Returns:
        List[List[Any]]: List of exon coordinates.
    """
    exon_data = re.findall(
        r"^(S\.\d+\s+\w+\s+[+-]\s+\d+\s+\d+|^\s*\d+\.\d+\s+\w+\s+[+-]\s+\d+\s+\d+)",
        results_text,
        flags=re.MULTILINE,
    )

    exon_list = []

    for exon in exon_data:
        if "Prom" not in exon and "PlyA" not in exon:
            exon_description = exon.split()
            gene_name = exon_description[0]
            exon_strand = exon_description[2]
            exon_start = exon_description[3]
            exon_end = exon_description[4]
            exon_list.append([gene_name, exon_strand, exon_start, exon_end])

    return exon_list


def extract_intron_list(exon_list: List[List[str]]) -> List[List[str]]:
    """
    Extracts intron coordinates from the list of exon coordinates.

    Args:
        exon_list (List[List[Any]]): List of exon coordinates.

    Returns:
        List[List[Any]]: List of intron coordinates.
    """
    intron_list = []
    for idx in range(len(exon_list) - 1):

        # take gene names for the pair of exons
        curr_exon_gene_name = exon_list[idx][0].split(".")[0]
        next_exon_gene_name = exon_list[idx + 1][0].split(".")[0]

        # take strands for the pair of exons
        curr_exon_strand = exon_list[idx][1]
        next_exon_strand = exon_list[idx + 1][1]

        # take exons in one gene on the same strand
        if (curr_exon_gene_name == next_exon_gene_name) and (
            curr_exon_strand == next_exon_strand
        ):
            # take full gene numbers of the exons pair and create intron name from it
            curr_exon_num = exon_list[idx][0]
            next_exon_num = exon_list[idx + 1][0]
            intron_gene_num = f"{curr_exon_num} - {next_exon_num}"

            # take same strand for intron
            intron_strand = curr_exon_strand

            # define intron start and end coordinates, based on strand direction
            if curr_exon_strand == "+":

                curr_exon_end = exon_list[idx][3]
                intron_start = int(curr_exon_end) + 1

                next_exon_start = exon_list[idx + 1][2]
                intron_end = int(next_exon_start) - 1

            else:

                next_exon_end = exon_list[idx + 1][3]
                intron_start = int(next_exon_end) + 1

                curr_exon_start = exon_list[idx][2]
                intron_end = int(curr_exon_start) - 1

            # collect output
            intron_list.append(
                [intron_gene_num, intron_strand, intron_start, intron_end]
            )

    return intron_list


def run_genscan(
    sequence: Optional[str] = None,
    sequence_file: Optional[str] = None,
    organism: str = "Vertebrate",
    exon_cutoff: float = 0.00,
    sequence_name: str = "",
) -> GenscanOutput:
    """
    Runs Genscan analysis.

    Args:
        sequence (Optional[str]): DNA sequence.
        sequence_file (Optional[str]): Path to a file containing DNA sequence.
        organism (str): Organism type for Genscan analysis.
        exon_cutoff (float): Exon cutoff value.
        sequence_name (str): Name of the DNA sequence.

    Returns:
        GenscanOutput: Object containing Genscan analysis results.
    """
    if sequence_file:
        sequence = ""
        if len(SeqIO.parse(sequence_file, "fasta")) > 1:
            print(
                "Your fasta file contains several records. Planning take first one for analysis."
            )

        for record in SeqIO.parse(sequence_file, "fasta")[:1]:
            sequence += str(record.seq)

    params = {
        "-s": sequence,
        "-o": organism,
        "-e": exon_cutoff,
        "-n": sequence_name,
        "-p": "Predicted peptides only",
    }

    response = requests.post(
        "http://hollywood.mit.edu/cgi-bin/genscanw_py.cgi", data=params
    )

    status = response.status_code

    if status != 200:
        raise Exception(
            f"Error: Unable to connect to Genscan service. Status code: {response.status_code}"
        )

    html_content = response.content
    soup = BeautifulSoup(html_content, "lxml")
    results_text = soup.find("pre").text

    cds_list = extract_cds_list(results_text)
    exon_list = extract_exon_list(results_text)
    intron_list = extract_intron_list(exon_list)

    return GenscanOutput(status, cds_list, intron_list, exon_list)
