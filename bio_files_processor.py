import os


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None) -> None:
    """
    Creates FASTA file with oneline sequences based on given FASTA file with multiline sequences in the same directory.
    :param input_fasta: path to the multiline FASTA file (str)
    :param output_fasta: name of output oneline FASTA file (str)
    :return: None
    This function creates a FASTA file in current directory.
    The new file will be named '<oneline_result_input_filename>.fasta' where <input_filename>
    is derived from the input file name by default or can be specified using output_fasta argument.
    """
    if output_fasta is None:
        output_fasta = f'oneline_result_{os.path.basename(input_fasta)}'
    if not output_fasta.endswith('.fasta'):
        output_fasta += '.fasta'

    path_to_out_dir = os.path.dirname(input_fasta)

    with open(input_fasta, 'r') as fin, open(os.path.join(path_to_out_dir, output_fasta)) as fout:
        seq = ''
        for line in fin:
            if line.startswith('>'):
                if seq:
                    fout.write(f'{seq}\n')
                    seq = ''
                fout.write(line)
            else:
                seq = f'{seq}\n'
        fout.write(seq)
