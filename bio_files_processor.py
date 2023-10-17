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


def select_genes_from_gbk_to_fasta(input_gbk: str, genes: list, n_before: int = 1, n_after: int = 1,
                                   output_fasta: str = None) -> None:
    """
    Creates FASTA file with neighbour CDSs to given genes from GBK file in fasta_selected_from_gbk directory.
    :param input_gbk: path to GBK file (str)
    :param genes: genes of interest that are used for searching neighbour CDSs (list)
    :param n_before: number of neighbor CDSs before gene of interest (int)
    :param n_after: number of neighbor CDSs after gene of interest (int)
    :param output_fasta: name of the output fasta file (str)
    :return: None
    """
    cds_coord_list = []
    coord_cds_dict = {}
    genes_names_dict = {}
    translation_seq = []
    record_translation = False
    gene_name = ''
    coord = ''
    with open(input_gbk, mode='r') as gbk:
        for line in gbk:
            if line.strip().startswith('CDS '):
                if coord:
                    coord_cds_dict[coord] = (gene_name, "".join(translation_seq))
                    translation_seq = []
                coord = line.split()[1]
                cds_coord_list.append(coord)
                record_translation = False
            elif '/gene' in line:
                gene_name = line.split('"')[1]
                genes_names_dict[gene_name] = coord
            if record_translation:
                translation_seq.append(line.strip().strip('"'))
            elif '/translation' in line:
                record_translation = True
                translation_seq.append(line.strip().split('"')[1])
            elif line.strip().startswith('ORIGIN'):
                record_translation = False
                coord_cds_dict[coord] = (gene_name, "".join(translation_seq))

        cds_of_interest = []

        for gene in genes:
            gene_coord = genes_names_dict[gene]
            gene_position = cds_coord_list.index(gene_coord)
            if len(cds_coord_list[:gene_position]) < n_before:
                raise ValueError("Too many neighbours CDSs before gene of interest are requested. Change number!")
            if len(cds_coord_list[gene_position:]) < n_after:
                raise ValueError("Too many neighbours CDSs after gene of interest are requested. Change number!")
            cds_of_interest = cds_coord_list[gene_position - n_before:gene_position + n_after + 1]

        os.makedirs('fasta_selected_from_gbk', exist_ok=True)

        if output_fasta is None:
            output_fasta = f"CDS_selected_from_gbk_{os.path.basename(input_gbk).split('.')[0]}"
        if not output_fasta.endswith('.fasta'):
            output_fasta = f'{output_fasta}.fasta'

        with open(os.path.join('fasta_selected_from_gbk', output_fasta), mode='w') as fasta:
            for cds in cds_of_interest:
                fasta.write(f'>{cds} gene:{coord_cds_dict[cds][0]}\n{coord_cds_dict[cds][1]}\n')
        print('FASTA file is ready.')
