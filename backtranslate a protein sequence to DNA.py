import itertools

def get_hexamita_codon_table():
    
    return {
        'F': ['TTT', 'TTC'],
        'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
        'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
        'Y': ['TAT', 'TAC'],
        'Q': ['TAA', 'TAG', 'CAA', 'CAG'],
        'C': ['TGT', 'TGC'],
        '*': ['TGA'],
        'W': ['TGG'],
        'P': ['CCT', 'CCC', 'CCA', 'CCG'],
        'H': ['CAT', 'CAC'],
        'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
        'I': ['ATT', 'ATC', 'ATA'],
        'M': ['ATG'],
        'T': ['ACT', 'ACC', 'ACA', 'ACG'],
        'N': ['AAT', 'AAC'],
        'K': ['AAA', 'AAG'],
        'V': ['GTT', 'GTC', 'GTA', 'GTG'],
        'A': ['GCT', 'GCC', 'GCA', 'GCG'],
        'D': ['GAT', 'GAC'],
        'E': ['GAA', 'GAG'],
        'G': ['GGT', 'GGC', 'GGA', 'GGG']
    }

def backtranslate_protein(protein_sequence):
    """
    Backtranslates a protein sequence to all possible DNA sequences
    based on the Hexamita Nuclear codon table.
    
    Args:
        protein_sequence (str): The protein sequence to backtranslate.
        
    Returns:
        list: A list of all possible DNA sequences.
    """
    codon_table = get_hexamita_codon_table()
    
    # Get a list of codon lists for each amino acid in the sequence
    codon_lists = [codon_table[aa.upper()] for aa in protein_sequence]
    
    # Use itertools.product to get all combinations of codons
    all_dna_combinations = list(itertools.product(*codon_lists))
    
    # Join the codons to form complete DNA sequences
    dna_sequences = [''.join(seq) for seq in all_dna_combinations]
    
    return dna_sequences

def analyze_codon_usage(protein_sequence):
    """
    Analyzes and reports the types of DNA states (codons) for each
    amino acid in a protein sequence.
    
    Args:
        protein_sequence (str): The protein sequence to analyze.
    
    Returns:
        None: Prints the analysis directly.
    """
    codon_table = get_hexamita_codon_table()
    print("Codon Usage Analysis for Hexamita Nuclear:")
    print("------------------------------------------")
    for aa in protein_sequence.upper():
        if aa in codon_table:
            codons = codon_table[aa]
            print(f"Amino Acid '{aa}' is encoded by {len(codons)} codon(s): {', '.join(codons)}")
        else:
            print(f"Amino Acid '{aa}' not found in the Hexamita Nuclear codon table.")


    print(seq)
