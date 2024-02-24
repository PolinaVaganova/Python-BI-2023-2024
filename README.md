# biopyseq-utils 
This repository contains an open-source library for processing protein, nucleic sequences and filter fastaq files. Our tool can process multiple proteins and peptides sequences, calculate physical features, find specific sites. For nucleic acids it can transcribe, find complement sequences, reverse, and many other things. Also it supports FASTAQ-format and can be used to filter DNA sequences based on their length, GC-content and sequence quality. 


## Installation

To use this toolbox one need to clone repository

```shell
git clone https://github.com/PolinaVaganova/pybioseq-utils/
cd pybioseq-utils
```
### System requirements:

Key packages and programs:
- [Python](https://www.python.org/downloads/) (version >= 3.9)


## biopyseq_utils.py

This script contains function which can:

- proccess nucleic seqs - **run_nucleic_seq_processing**

- proccess protein seqs - **run_protein_seq_processing**

- filter fastq file - **run_fastaq_filtering**

### Usage

```python
# import main functions
from biopyseq_utils import run_nucleic_seq_processing, run_protein_seq_processing, run_fastaq_filtering
```      
This section contains detailed description of biopyseq-utils functions.

### run_nucleic_seq_processing(\*args)

Apply one of the operations described below to any number of nucleic sequences.

**Parameters:**
**\*args**:
- **sequences**: *str*
input coma-separated dna or rna sequences in 1-letter code it can contain upper case and lower case letters
- **operation** : *str*
specify procedure you want to apply (see list below). Always must be the last argument!

**Returns**:
- **operation_result**: *str* or *list*
result of function work in list or str format (dependent on number of input sequences)

**Operations list**

- `transcribe` - transcribes given DNA sequences to RNA,

- `reverse_transcribe` - transcribe given RNA sequence to DNA,

- `reverse` - gives the reversed copy of a sequence

- `complement` - gives complement DNA or RNA sequence.

- `reverse_complement` - gives reversed complement DNA or RNA sequence. 

- `count_nucleotides` - gives percentages of each type of nucleotide in DNA or RNA sequence

- `make_triplets` - split your dna or rna seq into triplets

- `is_dna_or_rna` - check if sequence true nucleic sequence or not

**Example**
```python
dna_seq= 'ATCG'
rna_seq = 'AUUUGGC'

run_dna_rna_tools(dna_seq, 'transcribe')
run_dna_rna_tools(rna_seq, 'complement')
run_dna_rna_tools(dna_seq, rna_seq, 'count_nucleotides')
```

### run_protein_seq_processing(\*args)

Apply one of the operations described below to any number of protein sequences. 

**Parameters:**
**\*args**:
- **sequences**: *str*
input coma-separated protein sequences in 1-letter or 3-letter code
whitespaces may by in any place, because the will be deleted and replaced in output
sequences can contain upper case and lower case letters

- **current_encoding**: *str*
sequence encoding type in `str` format (`one` - one-letter encoding or `three` - three-letter encoding) 

- **add_arg**: *str*
optional parameter, required for certain operations (specify target protein site or residue of interest)
only 1-letter code is acceptable for residues in this argument!
argument can contain upper case and lower case letters

- **procedure** : *str*
specify procedure you want to apply (see list below). Always must be the last argument!

**Returns**:
- **operation_result**: *str*, *list*, *dict* or *float*
the result of procedure

**Note!**
- Arguments always must be in the strict order: sequences, sequence encoding type, additional argument, operation!

**Operations list**
- `change_residues_encoding` - transfer sequence in opposite from current tencoding

- `is_protein` - check if sequence is protein or not by identify invalid seq elements

- `get_seq_characteristic` - count entry of each residue type in your seq

- `find_residue_position` - find all positions of certain residue in your seq

- `find_site` - find if seq contains certain site and get positions of its site

- `calculate_protein_mass` - get mass of residues in your seq in Da

- `calculate_average_hydrophobicity` - get hydrophobicity index for protein seq as sum of index for each residue in your seq divided by its length

- `get_mrna` - get encoding mRNA nucleotides for your seq

- `calculate_isoelectric_point` - find isoelectrinc point as sum of known pI for residues in your seq

- `analyze_secondary_structure` - calculate the percentage of amino acids which responsible for the three main
    types of protein secondary structure: beta-turn, beta-sheet and alpha-helix

**Example**

```python
seq1 = 'Trp ala GLN phe'
seq2 = 'DAWLRHIL'
seq3 = 'YVWTD'

run_protein_analysis(seq1, 'three',  'AQF' 'find_site')
run_protein_analysis(seq2, 'one', 'get_mrna')
run_protein_analysis(seq2, seq3, 'one', 'D', 'find_residue_position')

```

### run_filter_fastaq(input_path, output_file, gc_bounds=(0, 100), length_bounds=(0, 2 ** 32), quality_threshold=0)

Filter any number of fasta files based on given on GC-content, length, and sequencing quality (phred33) from a FASTQ file.
Result will be saved in the `fastq_filtrator_results` directory within the same location as the input file.

**Parameters:**

**input_path**: *str*
- Path to the input FASTQ file

**output_filename** 
- Output FASTQ file name, which will be saved in `fastq_filtrator_results`. 
- If not provided, the output file will be named same as input FASTAQ file.

**gc_bounds**: *tuple*, *int* or *float*
- lower and upper boundaries of GC-content in percent orginised in *tuple*
- by default = (0, 100)
- if only one number given it will be considered the upper boundary with the lower boundary = 0

**length_bounds**: *tuple*, *int* or *float*
- lower and upper boundaries of sequence length orginised in *tuple*
- by default = (0, 2 ** 32)
- if only one number given it will be considered the upper boundary with the lover boundary = 0

**quality_threshold**: *int* or *float*
- lower boundary for threshold average Q-score
- by default=0

**Returns**:
- *None*
This function does not return anything. It saves the filtered FASTQ sequences
in the specified output file in `fastq_filtrator_results` directory.
    
**Example**
```python
run_filter_fastaq('/path/to/input/dir/input.fatq', output_filename='output.fastq', gc_bounds=(35, 82), length_bounds=(4, 8), quality_threshold=22)
# This will create a file named output.fastq in the fastq_filtrator_results directory with filtered sequences.
```

## bio_files_processor.py

### Usage

```python
# import main function
from bio_files_processor import convert_multiline_fasta_to_oneline, select_genes_from_gbk_to_fasta, change_fasta_start_pos, parse_blast_output

```
### convert_multiline_fasta_to_oneline(input_fasta, output_fasta) 

Creates FASTA file with oneline sequences based on given FASTA file with multiline sequences in the same directory.

**Parameters:**

**input_fasta**: *str*
- Path to the input multiline FASTA file

**output_fasta** : *str*
- Name of the output oneline FASTA file 
- If not provided, the output file will be named '<oneline_result_input_filename>.fasta'
    where <input_filename> is derived from the input file name.input FASTA  file.

**Returns**:
- *None*
This function creates the output FASTA file in current directory.
    
**Example**
```python
convert_multiline_fasta_to_oneline('/path/to/input/dir/input.fatq', output_fasta='output.fasta')
# This will create a file named output.fasta in the current directory.
```

### select_genes_from_gbk_to_fasta(input_gbk, genes, n_before, n_after, output_fasta)  

Extract neighbor CDSs to specified genes from a GBK file.
Generate a FASTA file containing the sequences in the `fasta_selected_from_gbk` directory within the same location as the input file.

**Parameters:**

**input_gbk**: *str*
- Path to the input GBK file

**genes** : *str*
- Genes of interest to be used for neighbor CDSs search 

**n_before** : *int*
- Number of neighbor CDSs before the gene of interest
- by default = 1

**n_after** : *int*
- Number of neighbor CDSs after the gene of interest
- by default = 1

**output_fasta**: *str*
- Name of the output FASTA file. If not provided, the output file will be named '<CDS_selected_from_gbk_input_filename>.fasta'
    where <input_filename> is derived from the input file name..

**Returns**:
- *None*
This function creates a FASTA file in `fasta_selected_from_gbk` directory.
    
**Example**
```python
select_genes_from_gbk_to_fasta('/path/to/input/dir/input.gbk', 'gene1', 'gene2', n_before=2, n_after=2, output_fasta='output.fasta')
# This will create a FASTA file named output.fasta in fasta_selected_from_gbk directory containing sequences of neighbor CDSs to 'gene1' and 'gene2' from the input.gbk file.
```

### change_fasta_start_pos(input_fasta, shift, output_fasta)  

Shift the starting position of sequence in the input FASTA file.

**Parameters:**

**input_fasta**: *str*
- Path to the input FASTA file

**shift** : *int*
- Number of positions to shift the start nucleotide of sequence

**output_fasta**: *str*
- Name of the output FASTA file. If not provided, the output file will be named 'shifted_by_<shift>_nucleotide_<input_filename>'
    where <input_filename> is derived from the input file name and <shift> from the corresponding param.

**Returns**:
- *None*
This function creates a FASTA file in current directory.

**Example**
```python
change_fasta_start_pos('/path/to/input/dir/input.fasta', 2, output_fasta='output.fasta')

# if input file contains `ATCCGT` sequence, the output.fasta will contain 'CCGTAT'

change_fasta_start_pos('/path/to/input/dir/input.fasta', 4, output_fasta='output.fasta')

# if input file contains `ATCCGT` sequence, the output.fasta will contain 'GTATCC'

```


### parse_blast_output(input_file, output_file)

Extract the descriptions (gene names) of the best BLAST results from the blast results file.
Result bill be saved in the `best_blast_results` directory within the same location as the input file.

**Parameters:**

**input_fasta**: *str*
- Path to the input FASTA file

**output_fasta**: *str*
- Name of the output txt file. If not provided, the output file will be named 'best_<input_filename>.txt', where
    <input_filename> is derived from the input file name.

**Returns**:
- *None*
This function creates a txt file in current directory.

**Example**
```python
parse_blast_output('/path/to/input/dir/blast_results.txt', output_file='output_results')

# This will create a file named output_results.txt with the descriptions of the best BLAST results.
```

Contribution
-----
Issues and PRs are welcome