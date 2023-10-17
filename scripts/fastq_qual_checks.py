from typing import Union


def gc_filtering(seq: str, gc_bounds: Union[tuple, float, int]) -> bool:
    """
    Filtering seq by GC content in percentages.
    :param seq: dna sequence
    :param gc_bounds: cutoff for GC content in percents. You can specify lower and upper limits (tuple with floats)
    or just upper limit (then pass float). Default=(0,100)
    :return: is seq in allowed zone (bool)
    """
    gc_percent = (seq.count('C') + seq.count('G')) / len(seq) * 100
    if isinstance(gc_bounds, tuple):
        return gc_bounds[0] <= gc_percent <= gc_bounds[1]
    return gc_percent <= gc_bounds


def length_filtering(seq: str, length_bounds: Union[tuple, float, int]) -> bool:
    """
    Filter seq by its length in nucleic bases.
    :param seq: dna sequence
    :param length_bounds: Lower and upper limits (tuple with floats)
    or just upper limit (float)
    :return: is seq in allowed zone (bool)
    """
    if isinstance(length_bounds, tuple):
        return length_bounds[0] <= len(seq) <= length_bounds[1]
    return len(seq) <= length_bounds


def quality_filtering(quality_seq: str, quality_threshold: Union[int, float]) -> bool:
    """
    Filter seq by its quality in phred33 scale.
    Reads with average score lower than cutoff will be dropped.
    :param quality_seq: quality of the sequence
    :param quality_threshold: cutoff for seq quality in phred33 scale.
    Reads with average score lower than this cutoff will be dropped
    :return: is seq passed cutoff (bool)
    """
    scores_list = []
    for quality_symbol in quality_seq:
        scores_list.append(ord(quality_symbol) - 33)
    return sum(scores_list) / len(scores_list) >= quality_threshold


def fastaq_to_dict(input_path: str) -> dict[str, tuple[str, str]]:
    """
    Parse input FASTQ file and write it content into dictionary
    :param input_path: path to the FASTQ file (str)
    :return: dict, where key is sequence name, values are seqs snd quality
    but the last key and value include input fastaq file name
    """
    pass


def dict_to_fastaq(seqs_dict: dict, output_filename: str = None) -> None:
    """
    Write dict with fastaq to the actual FASTQ file
    :param seqs_dict: dictionary with seqs and input FASTQ file name (dict)
    :param output_filename:  output FASTQ file name (str)
    :return: None
    This function creates a FASTQ file in 'fastq_filtrator_results' directory using the provided dictionary of sequences
    The new file will be named '<output_filename>.fastq' where <output_filename>
    is derived from the input file name by default or can be specified using output_filename argument.

    Example:
    --------
    If input file name is 'example.fastq' and output_filename is not specified the resulting FASTQ file will be
    'fastq_filtrator_results/example.fastq'.
    """
    pass
