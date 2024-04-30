import os

import pytest

from bio_files_processor import convert_multiline_fasta_to_oneline, OpenFasta
from pybioseq_utils import DNASequence, AminoAcidSequence, AminoAcidNotFoundError


class TestAminoAcidSequence:
    """
    Test class for AminoAcidSequence.
    """

    @pytest.fixture
    def input_data(self):
        return "AQEG"

    def test_check_seq(self):
        """
        Test to check the correctness of amino acid sequence.
        """
        with pytest.raises(AminoAcidNotFoundError):
            AminoAcidSequence("AQEGXXXX").check_seq_correctness()

    def test_count_aa(self, input_data):
        """
        Test to count amino acids in the sequence.
        """
        target = {"A": 1, "Q": 1, "E": 1, "G": 1}
        result = AminoAcidSequence(input_data).count_aa()
        assert target == result

    def test_str_aa(self, input_data):
        """
        Test string representation of AminoAcidSequence.
        """
        target = input_data
        result = str(AminoAcidSequence(input_data))
        assert target == result


class TestDNASequence:
    """
    Test class for DNASequence.
    """

    @pytest.fixture
    def input_data(self):
        return "ATAT"

    def test_dna_len(self, input_data):
        """
        Test to check the length of DNASequence.
        """
        target = 4
        result = len(DNASequence(input_data))
        assert target == result

    def test_dna_complement(self, input_data):
        """
        Test DNA complement.
        """
        target = "TATA"
        result = DNASequence(input_data).complement().seq
        assert target == result

    def test_dna_gc_content(self):
        """
        Test GC content of DNASequence.
        """
        target = 0
        result = DNASequence("ATAT").gc_content()
        assert target == result


class TestMultilineFastaProcessing:
    """
    Test class for convert_multiline_fasta_to_oneline function.
    """

    @pytest.fixture
    def input_data(self):
        return "example_data/example_fasta.fasta"

    @pytest.fixture
    def tmp_file(self):
        file_path = "tmp.fasta"
        yield file_path
        if os.path.exists(file_path):
            os.remove(file_path)

    def test_write_fasta_exists(self, input_data, tmp_file):
        """
        Test if the converted FASTA file exists.
        """
        convert_multiline_fasta_to_oneline(input_data, tmp_file)
        assert os.path.exists(tmp_file)

    def test_write_fasta_content(self, input_data, tmp_file):
        """
        Test the content of the converted FASTA file.
        """
        target = "ACGGCCATAGGACTTTGAAAGCACCGCATCCCGTCCGATCTGCGAAGTTAACCAAGATGCCGCCTGGTTAGTACCATGGTGGGGGACCACATGGGAATCCCTGGTGCTGTG"
        convert_multiline_fasta_to_oneline(input_data, tmp_file)
        with OpenFasta(tmp_file) as tmp:
            record = tmp.read_record()
            result = record.seq
        assert target == result
