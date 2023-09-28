def check_protein_seq(seq: str) -> str:
    
    """
    Checks whether a sequence is written using 1-letter amino acid code
    
    Arguments:
    -seq (str) input protein sequence
    Return:
    - str, 'single_letter_prot_seq' otherwise 'Invalid Input' error is raised
    
    """
    unique_chars = set(seq)
    single_letter = set('GALMFWKQESPVICYHRNDTgalmfwkqespvicyhrndt')

    if unique_chars <= single_letter:
        seq = 'single_letter_prot_seq'

    else:
        raise ValueError("Invalid Input")
    return seq
