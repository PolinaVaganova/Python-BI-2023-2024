from typing import Union, Sequence, List
import pybioseq_utils.fastaq as fastaqtutil
import pybioseq_utils.nucleic as nuclutil
import pybioseq_utils.protein as protutil


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
    function_names = {'transcribe': nuclutil.transcribe,
                      'reverse': nuclutil.reverse,
                      'complement': nuclutil.complement,
                      'reverse_complement': nuclutil.reverse_complement,
                      'count_nucleotides': nuclutil.count_nucleotides,
                      'make_triplets': nuclutil.make_triplets,
                      'reverse_transcribe': nuclutil.reverse_transcribe}
    procedure = args[-1]
    processed_result = []

    for ind, seq in enumerate(args[:-1]):
        if not nuclutil.is_dna_or_rna(seq):
            print(f'Sequence number {ind + 1} is not available for operations! Skip it.')
            continue
        processed_result.append(function_names[procedure](seq))
    if len(processed_result) == 1:
        return processed_result[0]
    return processed_result


# main function for nucleic seqs processing
def run_protein_seq_processing(*args: str) -> Union[List[str], str, float, List[float], dict, List[dict]]:
    """
    Specify and launch operation with proteins sequences. Parameters must be passed exactly in order given below

    :param args:
    - seq (str): amino acids sequences for analysis in 1-letter or 3-letter code (all encodings must be the same)
    Any number of sequences in any cases acceptable.
    Residues names may be separated by whitespaces (every 3 or 1 letter, depend on encoding).
    - additional arg (str): necessary parameter for some functions (for example, specify target protein site)
    - current_encoding (str): specify current encoding of your sequences as 'one' or 'three' corresponding to the length
    of residues names
    - operation name (str): specify procedure you want to apply

    :return: the result of procedure in list, str or float format
    """
    # first value is function name, second is real function, third is number of function additional arguments
    function_names = {'change_residues_encoding': [protutil.change_residues_encoding, 1],
                      'is_protein': [protutil.is_protein, 1],
                      'get_seq_characteristic': [protutil.get_seq_characteristic, 1],
                      'find_residue_position': [protutil.find_residue_position, 2],
                      'find_site': [protutil.find_site, 2],
                      'calculate_protein_mass': [protutil.calculate_protein_mass, 1],
                      'calculate_average_hydrophobicity': [protutil.calculate_average_hydrophobicity, 1],
                      'get_mrna': [protutil.get_mrna, 1],
                      'calculate_isoelectric_point': [protutil.calculate_isoelectric_point, 1],
                      'analyze_secondary_structure': [protutil.analyze_secondary_structure, 1]}

    # parse arguments
    procedure = args[-1]
    if function_names[procedure][1] == 1:
        current_encoding = args[-2]
    else:
        current_encoding = args[-3]
    seqs = [seq for seq in args[:-1 - (function_names[procedure][1])]]
    # create list for output
    processed_result = []
    # prepare sequence and launch desired operation
    for idx, seq in enumerate(seqs):
        if not protutil.is_protein(seq, current_encoding=current_encoding):
            print(f'Sequence number {idx + 1} is not available for operations! Skip it.')
            processed_result.append(f'Sequence number {idx + 1} is not available for operations! Skip it.')
            continue
        if current_encoding != 'one':
            seq = protutil.change_residues_encoding(seq, current_encoding=current_encoding)
        if function_names[procedure][1] == 1:
            if procedure == 'change_residues_encoding' or procedure == 'is_protein':
                processed_result.append(function_names[procedure][0](seq, current_encoding=current_encoding))
            else:
                processed_result.append(function_names[procedure][0](seq))
        elif function_names[procedure][1] == 2:
            add_arg = args[-2]
            processed_result.append(function_names[procedure][0](seq, add_arg))
    if len(processed_result) == 1:
        return processed_result[0]
    return processed_result


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
        gc_result = fastaqtutil.gc_filtering(seqs[seq_name][0], gc_bounds)
        length_result = fastaqtutil.length_filtering(seqs[seq_name][0], length_bounds)
        quality_result = fastaqtutil.quality_filtering(seqs[seq_name][1], quality_threshold)
        if gc_result and length_result and quality_result:
            filtered_seqs[seq_name] = seqs[seq_name]
    return filtered_seqs
