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

def molecular_weight(seq: str) -> int:
    
    """
    Calculates molecular weight of a protein
    
    Arguments:
    - seq (str) 1-letter coded protein sequence
    Return:
    - int, molecular weight (g/mol) rounded to integer   
        
    """
    if check_protein_seq(seq) == 'single_letter_prot_seq':
        list_input_seq = list(seq)
        aa_weight_dict = {'G':75, 'g':75, 'A':89, 'a':89, 'R':174, 'r':174, 'N':132, 'n':132, 
                          'D':133, 'd':133, 'C':121, 'c':133, 'E':147, 'e':147, 'Q':146, 'q':146,
                         'H':155, 'h':155, 'I':131, 'i':131, 'L':131, 'l':131, 'K':146, 'k':146, 
                         'M':149, 'm':149, 'F':165, 'f':165, 'P':115, 'p':115, 'S':105, 's':105, 
                         'T':119, 't':119, 'W':204, 'w':204, 'Y':181, 'y':181, 'V':117, 'v':117}
        water_mw = 18
        for aa in list_input_seq:
            total_mw = sum(aa_weight_dict[a] for a in list_input_seq)
            mw_water_removed = (total_mw - (water_mw * (len(list_input_seq)-1)))
        return mw_water_removed
