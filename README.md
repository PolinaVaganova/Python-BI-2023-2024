# biopyseq-utils 
This repository contains an open-source library for processing protein, nucleic sequences and filter fastaq files. Our tool can process multiple proteins and peptides sequences, calculate physical features, find specific sites. For nucleic acids it can transcribe, find complement sequences, reverse, and many other things. Also ir supports FASTAQ-format and can be used to filter DNA sequences based on their length, GC-content and sequence quality. 


## Installation

To use this toolbox one need to clone repository

```shell
git clone https://github.com/PolinaVaganova/pybioseq-utils/
cd pybioseq-utils
```
### System requirements:

Key packages and programs:
- [Python](https://www.python.org/downloads/) (version >= 3.9)

## Usage

```python
# import main functions
from biopyseq_utils import run_nucleic_seq_processing, run_protein_seq_processing, run_fastaq_filtering
```      
This section contains description of biopyseq-utils functions.

### run_nucleic_seq_processing(\*args)

Apply one of the operations described below to any number of nucleic sequences.

**Parameters:**
**\*args**:
- **sequences**: *str*
input coma-separated dna or rna sequences in 1-letter code it can contain upper case and lower case letters
- **operation** : *str*
specify procedure you want to apply (see list below). Always must be the last argument!

**Returns**:
- **operation_result**: str or list
result of function work in list or str format (dependent on number of input sequences)

**Operations list**

- `transcribe` - transcribes given DNA sequences to RNA,
	*Returns*
	`rna transcript: list[str] or str` 

- `reverse_transcribe` - transcribe given RNA sequence to DNA,
	*Returns*
	`dna seq: list[str] or str`

- `reverse` - gives the reversed copy of a sequence
	*Returns*
	`reverse seq: list[str] or str` 

- `complement` - gives complement DNA or RNA sequence.
	*Returns*
	`complement seq: list[str] or str` 

- `reverse_complement` - gives reversed complement DNA or RNA sequence. 
	*Returns*
	`result: list[str] or str` 

- `count_nucleotides` - gives percentages of each type of nucleotide in DNA or RNA sequence
	*Returns*
	`result: list[str] or str`

- `make_triplets` - split your dna or rna seq into triplets
	*Returns*
	`result: list[str] or str`

- `is_dna_or_rna` - check if sequence true nucleic sequence or not
	*Returns*
	`result: bool`

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
	*Returns*
	`seq in another code: list[str] or str` 

- `is_protein` - check if sequence is protein or not by identify invalid seq elements
	*Returns*:
	`result: list[bool] or bool`

- `get_seq_characteristic` - count entry of each residue type in your seq
	*Returns*
	`description of amino acid composition: dict`

- `find_residue_position` - find all positions of certain residue in your seq
	*Returns*
	`result: list[str] or str`

- `find_site` - find if seq contains certain site and get positions of its site
	*Returns*
	`result: list[str] or str`

- `calculate_protein_mass` - get mass of residues in your seq in Da
	*Returns*
	`mass: float or list[float]`

- `calculate_average_hydrophobicity` - get hydrophobicity index for protein seq as sum of index for each residue in your seq divided by its length
	*Returns*
	`hydrophobicity: float or list[float]`

- `get_mrna` - get encoding mRNA nucleotides for your seq
	*Returns*
	`mrna: str or list[str]`

- `calculate_isoelectric_point` - find isoelectrinc point as sum of known pI for residues in your seq
	*Returns*
	`pI: float or list[float]`

- `analyze_secondary_structure` - calculate the percentage of amino acids which responsible for the three main
    types of protein secondary structure: beta-turn, beta-sheet and alpha-helix
    *Returns*
    `result: str or list[str]`


**Example**

```python
seq1 = 'Trp ala GLN phe'
seq2 = 'DAWLRHIL'
seq3 = 'YVWTD'

run_protein_analysis(seq1, 'three',  'AQF' 'find_site')
run_protein_analysis(seq2, 'one', 'get_mrna')
run_protein_analysis(seq2, seq3, 'one', 'D', 'find_residue_position')

```

### run_fastaq_filtering(seqs, gc_bounds=(0, 100), length_bounds=(0, 2\*\*32), quality_threshold=0)

Filter any number of fasta files based on given arguments. 

**Parameters:**
**seqs**: *dict*
- FASTAQ sequences organised in dictionary: key = sequence name, value = sequence and sequence quality in ASCII code (all as *str*)

**gc_bounds**: *tuple*, *int* or *float*
- lower and upper boundaries of GC-content in percent orginised in `tuple`
- by default = (0, 100)
- if only one number given it will be considered the upper boundary with the lower boundary = 0

**length_bounds**: *tuple*, *int* or *float*
- lower and upper boundaries of sequence length orginised in `tuple`
- by default = (0, 2 ** 32)
- if only one number given it will be considered the upper boundary with the lover boundary = 0

**quality_threshold**: *int* or *float*
- lower boundary for threshold average Q-score
- by default=0

**Returns**:
- **filtered seq**: *dict*

**Returns**:
- **operation_result**: *str*, *list*, *dict* or *float*
Sequences which passed all the given cutoffs. 

**Example**
```python
fastaq_files_dict = {'seq_name1': ('ATCGATGCATGC', 'jjG#HHFddd@'), 'seq_name2': ('GGGTCATTT', '!@jHHj')}

run_filter_fastaq(fastaq_files_dict, gc_bounds=(30, 80), length_bounds=(4, 8), quality_threshold=22)

```

Contribution
-----
Issues and PRs are welcome