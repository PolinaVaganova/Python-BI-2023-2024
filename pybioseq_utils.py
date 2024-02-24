import os
from typing import Union
from Bio import SeqIO, SeqUtils


# function for fastaq files filtering
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
