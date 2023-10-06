from typing import Union, Sequence, List


def gc_filtering(seq: str, gc_bounds: Union[tuple, float, int]) -> bool:
    """
    Filtering seq by GC content in percentages.
    :param seq: dna sequence
    :param gc_bounds: cutoff for GC content in percents. You can specify lower and upper limits (tuple with floats)
    or just upper limit (then pass float). Default=(0,100)
    :return: is seq in allowed zone (bool)
    """
    gc_percent = (seq.count('C') + seq.count('G')) / len(seq) * 100
    if type(gc_bounds) == tuple:
        return gc_percent < gc_bounds[0] or gc_percent > gc_bounds[1]
    else:
        return gc_percent < gc_bounds



def length_filtering(seq: str, length_bounds: tuple) -> bool:
    """
    Filter seq by its length in nucleic bases.
    :param seq: dna sequence
    :param length_bounds: Lower and upper limits (tuple with floats)
    or just upper limit (float)
    :return: is seq in allowed zone (bool)
    """
    pass


def quality_filtering(quality_seq: str, quality_threshold: int) -> bool:
    """
    Filter seq by its quality in phred33 scale.
    Reads with average score lower than cutoff will be dropped.
    :param quality_seq: quality of the sequence
    :param quality_threshold: cutoff for seq quality in phred33 scale.
    Reads with average score lower than this cutoff will be dropped
    :return: is seq passed cutoff (bool)
    """
    pass


print(gc_filtering('GCGCGCGC', 20))
