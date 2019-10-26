aa_scales = {

    # Charge
    'charge': {'K': 1, 'R': 1, 'D': -1, 'E': -1},
    'positive_charge': {'K': 1, 'R': 1},
    'negative_charge': {'D': -1, 'E': -1},
    
    # From BioPython Bio.SeqUtils.ProtParam
    'monoisotopic_protein_weights': {'A': 89.047678, 'C': 121.019749, 'D': 133.037508, 'E': 147.053158, 'F': 165.078979, 'G': 75.032028, 'H': 155.069477, 'I': 131.094629, 'K': 146.105528, 'L': 131.094629, 'M': 149.051049, 'N': 132.053492, 'O': 255.158292, 'P': 115.063329, 'Q': 146.069142, 'R': 174.111676, 'S': 105.042593, 'T': 119.058243, 'U': 168.964203, 'V': 117.078979, 'W': 204.089878, 'Y': 181.073893},
    'aromaticity': {'Y': 1, 'W': 1, 'F': 1},
    'gravy': {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y':-1.3, 'V': 4.2},
    'ss_helix': {'V': 1, 'I' : 1, 'Y': 1, 'F': 1, 'W': 1, 'L': 1},
    'ss_turn': {'N': 1, 'P': 1, 'G': 1, 'S': 1},
    'ss_sheet': {'E': 1, 'M': 1, 'A': 1, 'L': 1},
    
    # Kyte & Doolittle index of hydrophobicity
    'index_of_hydrophobicity': {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}, 
    
    # Hydrophilicity: 1 Hopp & Wood, Proc. Natl. Acad. Sci. U.S.A. 78: 3824-3828(1981)
    'hydrophilicity': {'A': -0.5, 'R': 3.0, 'N': 0.2, 'D': 3.0, 'C': -1.0, 'Q': 0.2, 'E': 3.0, 'G': 0.0, 'H': -0.5, 'I': -1.8, 'L': -1.8, 'K': 3.0, 'M': -1.3, 'F': -2.5, 'P': 0.0, 'S': 0.3, 'T': -0.4, 'W': -3.4, 'Y': -2.3, 'V': -1.5},
    
    # Normalized flexibility parameters (B-values) average (Vihinen et al., 1994)
    'flexibility': {'A': 0.984, 'C': 0.906, 'E': 1.094, 'D': 1.068, 'G': 1.031, 'F': 0.915, 'I': 0.927, 'H': 0.950, 'K': 1.102, 'M': 0.952, 'L': 0.935, 'N': 1.048, 'Q': 1.037, 'P': 1.049, 'S': 1.046, 'R': 1.008, 'T': 0.997, 'W': 0.904, 'V': 0.931, 'Y': 0.929}, 
    
    # Surface accessibility {"em"}: 1 Emini Surface fractional probability
    'surface_accessibility': {'A': 0.815, 'R': 1.475, 'N': 1.296, 'D': 1.283, 'C': 0.394, 'Q': 1.348, 'E': 1.445, 'G': 0.714, 'H': 1.180, 'I': 0.603, 'L': 0.603, 'K': 1.545, 'M': 0.714, 'F': 0.695, 'P': 1.236, 'S': 1.115, 'T': 1.184, 'W': 0.808, 'Y': 1.089, 'V': 0.606}, 
    
    # 2 Janin Interior to surface transfer energy scale
    'surface_transfer_energy': {'A': 0.28, 'R': -1.14, 'N': -0.55, 'D': -0.52, 'C': 0.97, 'Q': -0.69, 'E': -1.01, 'G': 0.43, 'H': -0.31, 'I': 0.60, 'L': 0.60, 'K': -1.62, 'M': 0.43, 'F': 0.46, 'P': -0.42, 'S': -0.19, 'T': -0.32, 'W': 0.29, 'Y': -0.15, 'V': 0.60}, 
    
    # Positive values indicate proteins that are likely to be ordered
    'ordered': {'A': 0.06, 'R': 0.180, 'N': 0.007, 'D': 0.192, 'C': 0.02, 'Q': 0.318, 'E': 0.736, 'G': 0.166, 'H': 0.303, 'I': -0.486, 'L': -0.326, 'K': 0.586, 'M': -0.397, 'F': -0.697, 'P': 0.987, 'S': 0.341, 'T': 0.059, 'W': -0.884, 'Y': -0.510, 'V': -0.121}, 

    # https: //github.com/ddofer/Protein-Descriptors/blob/master/src/csdsML/Descriptors.py
    'polarizability': {'A': 0.046, 'R': 0.291, 'N': 0.134, 'D': 0.105, 'C': 0.128, 'Q': 0.180, 'E': 0.151, 'G': 0.000, 'H': 0.230, 'I': 0.186, 'L': 0.186, 'K': 0.219, 'M': 0.221, 'F': 0.290, 'P': 0.131, 'S': 0.062, 'T': 0.108, 'W': 0.409, 'Y': 0.298, 'V': 0.140}, 
    'asa_in_tripeptide': {'A': 115, 'R': 225, 'N': 160, 'D': 150, 'C': 135, 'Q': 180, 'E': 190, 'G': 75, 'H': 195, 'I': 175, 'L': 170, 'K': 200, 'M': 185, 'F': 210, 'P': 145, 'S': 115, 'T': 140, 'W': 255, 'Y': 230, 'V': 155}, 
    'volume': {'A': 52.6, 'R': 109.1, 'N': 75.7, 'D': 68.4, 'C': 68.3, 'Q': 89.7, 'E': 84.7, 'G': 36.3, 'H': 91.9, 'I': 102.0, 'L': 102.0, 'K': 105.1, 'M': 97.7, 'F': 113.9, 'P': 73.6, 'S': 54.9, 'T': 71.2, 'W': 135.4, 'Y': 116.2, 'V': 85.1}, 
    'steric_param': {'A': 0.52, 'R': 0.68, 'N': 0.76, 'D': 0.76, 'C': 0.62, 'Q': 0.68, 'E': 0.68, 'G': 0.00, 'H': 0.70, 'I': 1.02, 'L': 0.98, 'K': 0.68, 'M': 0.78, 'F': 0.70, 'P': 0.36, 'S': 0.53, 'T': 0.50, 'W': 0.70, 'Y': 0.70, 'V': 0.76},
    'mutability': {'A': 100, 'R': 65, 'N': 134, 'D': 106, 'C': 20, 'Q': 93, 'E': 102, 'G': 49, 'H': 66, 'I': 96, 'L': 40, 'K': 56, 'M': 94, 'F': 41, 'P': 56, 'S': 120, 'T': 97, 'W': 18, 'Y': 41, 'V': 74},
}

ordered_aa_scale_names = ['gravy', 'ss_helix', 'ordered', 'mutability', 'steric_param', 'polarizability', 'flexibility', 'negative_charge', 'hydrophilicity', 'ss_turn', 'volume', 'surface_accessibility', 'charge', 'index_of_hydrophobicity', 'surface_transfer_energy', 'asa_in_tripeptide', 'ss_sheet', 'monoisotopic_protein_weights', 'aromaticity', 'positive_charge']