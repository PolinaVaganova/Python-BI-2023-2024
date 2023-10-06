from typing import Union, Sequence, List
import pybioseq_utils.fastaq as fastaqt
import pybioseq_utils.nucleic
import pybioseq_utils.protein


# main function for nucleic seqs processing
def run_nucleic_seq_processing(*args: str) -> Union[List[Sequence], Sequence]:
    """
    Specify and launch operation with dna or rna sequences

    :param args:
    - seq (str): dna or rna sequences for analysis (any number and case)
    - operation name (str): chosen procedure for analysis

    :return:
    - str: the result of procedure (case-sensitive)
    """
    pass


# main function for nucleic seqs processing
def run_protein_seq_processing(*args: str) -> Union[List[str], str, float, List[float]]:
    """
    Specify and launch operation with proteins sequences.

    :param args:
    - seq (str): amino acids sequences for analysis in 1-letter or 3-letter code
    Any number and cases. Whitespaces might be only between every residue name (every 3 or 1 letter, depend on encoding)
    - additional arg (str): necessary parameter for certain functions (for example, specify target protein site)
    - operation name (str): specify procedure you want to apply

    :return: the result of procedure in list, str or float format
    """
    pass


# main function for fastaq seqs filtering
def run_fastaq_filtering(seqs: dict, gc_bounds: Union[tuple, float, int] = (0, 100),
                         length_bounds: Union[tuple, float] = (0, 2 ** 32),
                         quality_threshold: int = 0) -> dict:
    """
    Launch filtering fastaq seq using 3 adjustable cutoffs. Allowed intervals include cutoffs values.

    :param seqs: dict with fastaq seqs, where key - seq name (str), value - seq, it's quality (tuple with str).
    :param gc_bounds: cutoff for GC content in percents. You can specify lower and upper limits (tuple with floats)
    or just upper limit (then pass float). Default = (0,100)
    :param length_bounds: cutoff fot length in nucleic bases. You can specify lower and upper limits (tuple with floats)
    or just upper limit (then pass float ot int). Default = (0,2**32)
    :param quality_threshold: cutoff for seq quality in phred33 scale. Default = 0.
    Reads with average score lower than this cutoff will be dropped.

    :return:
    """
    filtered_seqs = dict()
    for seq_name in seqs:
        gc_result = fastaqt.gc_filtering(seqs[seq_name][0], gc_bounds)
        length_result = fastaqt.length_filtering(seqs[seq_name][0], length_bounds)
        quality_result = fastaqt.quality_filtering(seqs[seq_name][1], quality_threshold)
        if gc_result and length_result and quality_result:
            filtered_seqs[seq_name] = seqs[seq_name]
    return filtered_seqs


