import os
from dataclasses import dataclass
from typing import List, Any


class OpenFasta:
    """
    Represents a context manager for reading FASTA files.

    Attributes:
        file_path (str): The path to the FASTA file.
        iterator (FastaIterator): Iterator object for iterating over records in the FASTA file.

    Methods:
        __enter__(): Opens the FASTA file and returns the OpenFasta object.
        __iter__(): Returns the iterator object for iterating over records in the FASTA file.
        read_record(): Reads a single record from the FASTA file.
        read_records(): Reads all records from the FASTA file.
        __exit__(): Closes the file when exiting the context manager.
    """

    def __init__(self, file_path: str):
        """
        Initializes an OpenFasta object.

        Args:
            file_path (str): The path to the FASTA file.
        """
        self.file_path = file_path
        self.iterator = None

    def __enter__(self) -> "OpenFasta":
        """
        Opens the FASTA file and returns the OpenFasta object.

        Returns:
            OpenFasta: The OpenFasta object.
        """
        self.file_handler = open(self.file_path, mode="r")
        return self

    def __iter__(self) -> "FastaIterator":
        """
        Returns the iterator object for iterating over records in the FASTA file.

        Returns:
            FastaIterator: The iterator object.
        """
        if self.iterator is None:
            self.iterator = FastaIterator(self.file_handler)
        return self.iterator

    def read_record(self) -> Any:
        """
        Reads a single record from the FASTA file.

        Returns:
            Any: A single record from the FASTA file.
        """
        return self.__iter__().__next__()

    def read_records(self) -> List[Any]:
        """
        Reads all records from the FASTA file.

        Returns:
            List[Any]: A list of all records from the FASTA file.
        """
        all_records = []
        for record in self.__iter__():
            all_records.append(record)
        return all_records

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Closes the file when exiting the context manager.

        Args:
            exc_type: Exception type.
            exc_val: Exception value.
            exc_tb: Exception traceback.
        """
        if self.file_handler:
            self.file_handler.close()


@dataclass
class FastaRecord:
    """
    Represents a single FASTA record.

    Attributes:
        id (str): The identifier of the record.
        description (str): The description of the record.
        seq (str): The sequence data of the record.

    Methods:
        __repr__(): Returns a string representation of the FASTA record.
    """

    id: str
    description: str
    seq: str

    def __repr__(self) -> str:
        """
        Returns a string representation of the FASTA record.

        Returns:
            str: A string representation of the FASTA record.
        """
        header = f">{self.id} {self.description}\n"
        return f"{header}{self.seq}"


class FastaIterator:
    """
    Iterator for iterating over records in a FASTA file.

    Attributes:
        _ids (List[str]): List of record identifiers.
        _seqs (List[str]): List of sequence data.
        _descriptions (List[str]): List of record descriptions.
        current_pos (int): Current position in the iterator.
        _open_fasta: File handler for the open FASTA file.
    """

    def __init__(self, open_fasta: "OpenFasta"):
        """
        Initializes a FastaIterator object.

        Args:
            open_fasta (OpenFasta): The OpenFasta object.
        """
        self._ids = []
        self._seqs = []
        self._descriptions = []

        self.current_pos = 0
        self._open_fasta = open_fasta

    def __iter__(self) -> "FastaIterator":
        """
        Returns the iterator object.

        Returns:
            FastaIterator: The iterator object.
        """
        return self

    def __next__(self) -> FastaRecord:
        """
        Iterates over records in the FASTA file.

        Returns:
            FastaRecord: The next FASTA record.

        Raises:
            StopIteration: If there are no more records to read.
        """
        if self.current_pos == len(self._seqs):
            seq = ""
            line = self._open_fasta.readline()
            while line:
                if line.startswith(">"):
                    id_and_description = line[1:].strip().split(" ", 1)
                    self._ids.append(id_and_description[0])
                    self._descriptions.append(id_and_description[1])
                    if seq:
                        self._seqs.append(seq)
                        seq = ""
                else:
                    seq += line.strip()
                line = self._open_fasta.readline()

            if seq:
                self._seqs.append(seq)
            if self.current_pos == len(self._seqs):
                raise StopIteration

        current_pos = self.current_pos
        self.current_pos += 1
        return FastaRecord(
            self._ids[current_pos],
            self._descriptions[current_pos],
            self._seqs[current_pos],
        )


def convert_multiline_fasta_to_oneline(
    input_fasta: str, output_fasta: str = None
) -> None:
    """
    Creates FASTA file with oneline sequences based on given FASTA file with multiline sequences.
    :param input_fasta: path to the multiline FASTA file (str)
    :param output_fasta: name of output oneline FASTA file (str)
    :return: None
    This function creates a FASTA file in current directory.
    If the output_file param is not specified the new file will be named '<oneline_result_input_filename>.fasta'
    where <input_filename> is derived from the input file name.
    """
    if output_fasta is None:
        output_fasta = f"oneline_result_{os.path.basename(input_fasta)}"
    if not output_fasta.endswith(".fasta"):
        output_fasta += ".fasta"

    with open(input_fasta, "r") as fin, open(output_fasta, "w") as fout:
        seq = ""
        for line in fin:
            if line.startswith(">"):
                if seq:
                    fout.write(f"{seq}\n")
                    seq = ""
                fout.write(line)
            else:
                seq = f"{seq}\n"
        fout.write(seq)


def select_genes_from_gbk_to_fasta(
    input_gbk: str,
    genes: list,
    n_before: int = 1,
    n_after: int = 1,
    output_fasta: str = None,
) -> None:
    """
    Creates FASTA file with neighbour CDSs to given genes from GBK file in fasta_selected_from_gbk directory.
    :param input_gbk: path to GBK file (str)
    :param genes: genes of interest that are used for searching neighbour CDSs (list)
    :param n_before: number of neighbor CDSs before gene of interest (int)
    :param n_after: number of neighbor CDSs after gene of interest (int)
    :param output_fasta: name of the output fasta file (str)
    :return: None
    This function creates a FASTA file in fasta_selected_from_gbk directory.
    If the output_file param is not specified the new file will be named '<CDS_selected_from_gbk_input_filename>.fasta'
    where <input_filename> is derived from the input file name.
    """
    cds_coord_list = []
    coord_cds_dict = {}
    genes_names_dict = {}
    translation_seq = []
    record_translation = False
    gene_name = ""
    coord = ""
    with open(input_gbk, mode="r") as gbk:
        for line in gbk:
            if line.strip().startswith("CDS "):
                if coord:
                    coord_cds_dict[coord] = (gene_name, "".join(translation_seq))
                    translation_seq = []
                coord = line.split()[1]
                cds_coord_list.append(coord)
                record_translation = False
            elif "/gene" in line:
                gene_name = line.split('"')[1]
                genes_names_dict[gene_name] = coord
            if record_translation:
                translation_seq.append(line.strip().strip('"'))
            elif "/translation" in line:
                record_translation = True
                translation_seq.append(line.strip().split('"')[1])
            elif line.strip().startswith("ORIGIN"):
                record_translation = False
                coord_cds_dict[coord] = (gene_name, "".join(translation_seq))

        cds_of_interest = []

        for gene in genes:
            gene_coord = genes_names_dict[gene]
            gene_position = cds_coord_list.index(gene_coord)
            if len(cds_coord_list[:gene_position]) < n_before:
                raise ValueError(
                    "Too many neighbours CDSs before gene of interest are requested. Change number!"
                )
            if len(cds_coord_list[gene_position:]) < n_after:
                raise ValueError(
                    "Too many neighbours CDSs after gene of interest are requested. Change number!"
                )
            cds_of_interest = cds_coord_list[
                gene_position - n_before : gene_position + n_after + 1
            ]

        os.makedirs("fasta_selected_from_gbk", exist_ok=True)

        if output_fasta is None:
            output_fasta = (
                f"CDS_selected_from_gbk_{os.path.basename(input_gbk).split('.')[0]}"
            )
        if not output_fasta.endswith(".fasta"):
            output_fasta = f"{output_fasta}.fasta"

        with open(
            os.path.join("fasta_selected_from_gbk", output_fasta), mode="w"
        ) as fasta:
            for cds in cds_of_interest:
                fasta.write(
                    f">{cds} gene:{coord_cds_dict[cds][0]}\n{coord_cds_dict[cds][1]}\n"
                )
        print("FASTA file is ready.")


def change_fasta_start_pos(input_fasta: str, shift: int, output_fasta=None) -> None:
    """
    Shift the starting position of sequence in the input FASTA file
    :param input_fasta: path to the input FASTA file (str)
    :param shift: number of positions to shift the start nucleotide of sequence (int)
    :param output_fasta: name of the output FASTA file (str)
    :return: None
    This function creates a FASTA file in current directory.
    If the output_file param is not specified the new file will be named
    '<shifted_by_shift_nucleotide_input_filename>.fasta' where <input_filename>
    is derived from the input file name.

    """
    if output_fasta is None:
        output_fasta = f"shifted_by_{shift}_nucleotide_{os.path.basename(input_fasta)}"
    if not output_fasta.endswith(".fasta"):
        output_fasta = f"{output_fasta}.fasta"
    with open(input_fasta, mode="r") as fin, open(output_fasta, mode="w") as fout:
        for line in fin:
            if line.startswith(">"):
                fout.write(f"{line}\n")
            else:
                fout.write(f"{line[shift:]}{line[:shift]}")


def parse_blast_output(input_file: str, output_file=None) -> None:
    """
    Write descriptions of the best blast results to the file.
    :param input_file: path to input blast results file (str)
    :param output_file: name for output file (str);
    This function creates a txt file in current directory.
    :return: None
    This function creates a txt file in current directory.
    If the output_file param is not specified the new file will be named 'best_<input_filename>.txt', where
    <input_filename> is derived from the input file name.

    """
    best_blast_results = []
    with open(input_file, mode="r") as fin:
        query = False
        description = False
        for line in fin:
            if line.startswith("Query #"):
                query = True
            elif query and line.startswith("Description  "):
                description = True
            elif query and description:
                current_result = line.split("    ")[0]
                best_blast_results.append(current_result)
                query = False
                description = False
    if output_file is None:
        output_file = f"best_{os.path.basename(input_file).split('.')[0]}.txt"
    if not output_file.endswith(".txt"):
        output_file = f"{output_file}.fasta"
    with open(output_file, mode="w") as fout:
        for description in best_blast_results:
            fout.write(f"{description}\n")
