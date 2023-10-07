from typing import List, Union

# 3-letter with corresponding 1-letter residues names
RESIDUES_NAMES = {'ALA': 'A',
                  'ARG': 'R',
                  'ASN': 'N',
                  'ASP': 'D',
                  'CYS': 'C',
                  'GLN': 'Q',
                  'GLU': 'E',
                  'GLY': 'G',
                  'HIS': 'H',
                  'ILE': 'I',
                  'LEU': 'L',
                  'LYS': 'K',
                  'MET': 'M',
                  'PHE': 'F',
                  'PRO': 'P',
                  'SER': 'S',
                  'THR': 'T',
                  'TRP': 'W',
                  'TYR': 'Y',
                  'VAL': 'V'
                  }

# first value is hydrophobicity index, second is pKa (pKa1, pKa2, pKa3 respectively), third is molecular mass in Da
RESIDUES_CHARACTERISTICS = {'A': [1.8, [2.34, 9.69, 0], 89],
                            'R': [-4.5, [2.17, 9.04, 12.48], 174],
                            'N': [-3.5, [2.02, 8.80, 0], 132],
                            'D': [-3.5, [1.88, 9.60, 3.65], 133],
                            'C': [2.5, [1.96, 10.28, 8.18], 121],
                            'Q': [-3.5, [2.17, 9.13, 0], 146],
                            'E': [-3.5, [2.19, 9.67, 4.25], 147],
                            'G': [-0.4, [2.34, 9.60, 0], 75],
                            'H': [-3.2, [1.82, 9.17, 6.00], 155],
                            'I': [4.5, [2.36, 9.60, 0], 131],
                            'L': [3.8, [2.36, 9.60, 0], 131],
                            'K': [-3.9, [2.18, 8.95, 10.53], 146],
                            'M': [1.9, [2.28, 9.21, 0], 149],
                            'F': [2.8, [1.83, 9.13, 0], 165],
                            'P': [-1.6, [1.99, 10.60, 0], 115],
                            'S': [-0.8, [2.21, 9.15, 0], 105],
                            'T': [-0.7, [2.09, 9.10, 0], 119],
                            'W': [-0.9, [2.83, 9.39, 0], 204],
                            'Y': [-1.3, [2.20, 9.11, 0], 181],
                            'V': [4.2, [2.32, 9.62, 0], 117]}

# amino acid with corresponding degenerate codon/codons
AMINO_ACID_TO_MRNA = {'A': 'GCN',
                      'R': '(CGN/AGR)',
                      'N': 'AAY',
                      'D': 'GAY',
                      'C': 'UGY',
                      'Q': 'CAR',
                      'E': 'GAR',
                      'G': 'GGN',
                      'H': 'CAY',
                      'I': 'AUH',
                      'L': '(CUN/UUR)',
                      'K': 'AAR',
                      'M': 'AUG',
                      'F': 'UUY',
                      'P': 'CCN',
                      'S': '(UCN/AGY)',
                      'T': 'ACN',
                      'W': 'UGG',
                      'Y': 'UAY',
                      'V': 'GUN'}


def change_residues_encoding(seq: str, current_encoding: str) -> str:
    """
    Transfer amino acids from 3-letter to 1-letter code and vice versa. By default, converts all seq into 1-letter
    format, even those already 1-letter. Case-insensitive.
    :param seq: protein seq (str)
    :param current_encoding: specify current encoding (str)
    :return: same protein seq in another encoding (str)
    """
    seq_in_target_encoding = []
    temp_seq = seq.replace(' ', '')
    if current_encoding == 'three':
        for residue_start_idx in range(0, len(temp_seq) - 2, 3):
            residue = temp_seq[residue_start_idx:residue_start_idx + 3]
            seq_in_target_encoding.append(RESIDUES_NAMES[residue.upper()])
    elif current_encoding == 'one':
        for residue in temp_seq:
            seq_in_target_encoding.append([target_residue for target_residue in RESIDUES_NAMES if
                                           RESIDUES_NAMES[target_residue] == residue.upper()][0])
            seq_in_target_encoding.append(' ')
    else:
        raise ValueError("Please, specify current encoding of your sequences as 'one' or 'three'")

    return ''.join(seq_in_target_encoding)


def is_protein(seq: str, current_encoding: str) -> bool:
    """
    Check if sequence is protein or not by identify invalid seq elements, which are not presented in dicts above.
    :param seq: protein seq (str)
    :param current_encoding: specify current encoding (str)
    :return: if seq is correct protein seq or not (bool)
    """
    temp_seq = seq.replace(' ', '')
    if current_encoding == 'one':
        for residue in temp_seq:
            if residue.upper() not in RESIDUES_NAMES.values():
                return False
        return True
    if current_encoding == 'three':
        for residue_start_idx in range(0, len(temp_seq) - 2, 3):
            residue = temp_seq[residue_start_idx:residue_start_idx + 3]
            if residue.upper() not in RESIDUES_NAMES.values():
                return False
        return True


def get_seq_characteristic(seq: str) -> dict:
    """
    Count entry of each residue type in your seq. Get description of amino acid composition in dict format.
    :param seq: protein seq in 1-letter encoding (str)
    :return: each residue type in seq in 3-letter code and its amount in current seq (dict)
    """
    res_count = {}
    for residue in set(seq):
        residue_entry = seq.count(residue)
        res_count[[tl_code for tl_code in RESIDUES_NAMES if RESIDUES_NAMES[tl_code] == residue][0]] = residue_entry
    return res_count


def find_residue_position(seq: str, res_of_interest: str) -> str:
    """
    Find all positions of certain residue in your seq
    :param seq: protein seq in 1-letter encoding (str)
    :param res_of_interest: specify the residue of interest in 1-letter encoding (str)
    :return: positions of specified residue in your seq (str)
    """
    res_of_interest = res_of_interest.upper()
    res_of_interest_positions = []
    for ind, res in enumerate(seq, 1):
        if res == res_of_interest:
            res_of_interest_positions.append(ind)
    return f"{res_of_interest} positions: {', '.join(map(str, res_of_interest_positions))}"


def find_site(seq: str, site: str) -> str:
    """
    Find if seq contains certain site and get positions of its site
    :param seq: protein seq in 1-letter encoding (str)
    :param site: specify site of interest in 1-letter encoding (str)
    :return: positions of residues for each certain site in seq (str)
    """
    site = site.upper()
    if not is_protein(site, current_encoding='one'):
        return f'Site {site} is not a protein!'
    if site in seq:
        site_full_position = []
        site_count = seq.count(site)
        site_start_position = [(coordinate + 1) for coordinate in range(len(seq)) if seq.startswith(site, coordinate)]
        site_end_position = [(coordinate + len(site) - 1) for coordinate in site_start_position]
        for idx in range(len(site_start_position)):
            site_full_position.append(f'{site_start_position[idx]}:{site_end_position[idx]}')
        return f"Site entry in sequence = {site_count}. Site residues can be found at positions: {', '.join(site_full_position)}"
    else:
        return f"{site} site is not in sequence!"


def calculate_protein_mass(seq: str) -> float:
    """
    Get mass of residues in your seq in Da
    :param seq: protein seq in 1-letter encoding (str)
    :return: mass in Da (float)
    """
    total_mass = 0
    for res in seq.upper():
        total_mass += RESIDUES_CHARACTERISTICS[res][2]
    return total_mass


def calculate_average_hydrophobicity(seq: str) -> float:
    """
    Get hydrophobicity index for protein seq as sum of index for each residue in your seq divided by its length
    :param seq: protein seq in 1-letter encoding (str)
    :return: average hydrophobicity (float)
    """
    sum_hydrophobicity_ind = 0
    for res in seq.upper():
        sum_hydrophobicity_ind += RESIDUES_CHARACTERISTICS[res][0]
    return sum_hydrophobicity_ind / len(seq)


def get_mrna(seq: str) -> str:
    """
    Get encoding mRNA nucleotides for your seq
    :param seq: protein seq in 1-letter encoding (str)
    :return: potential encoding mRNA sequence with multiple choice for some positions (str)
    """
    mrna_seq = str()
    for res in seq.upper():
        mrna_seq += AMINO_ACID_TO_MRNA[res]
    return mrna_seq


def calculate_isoelectric_point(seq: str) -> float:
    """
    Find isoelectrinc point as sum of known pI for residues in your seq
    :param seq: protein seq in 1-letter encoding (str)
    :return: isoelectric point (float)
    """
    sum_pka = 0
    pka_amount = 0
    for ind, res in enumerate(seq.upper(), 1):
        if ind == 1:
            sum_pka += RESIDUES_CHARACTERISTICS[res][1][1]
            pka_amount += 1
        elif RESIDUES_CHARACTERISTICS[res][1][2] != 0:
            sum_pka += RESIDUES_CHARACTERISTICS[res][1][2]
            pka_amount += 1
        elif ind == len(seq):
            sum_pka += RESIDUES_CHARACTERISTICS[res][1][0]
            pka_amount += 1
    pi = sum_pka / pka_amount
    return pi


def analyze_secondary_structure(seq: str) -> list[str]:
    """
    Calculate the percentage of amino acids found in the three main
    types of protein secondary structure: beta-turn, beta-sheet and alpha-helix
    :param seq: protein seq in 1-letter encoding (str)
    :return: percentage of amino acids belonging to three types of secondary structure (list[str])
    """
