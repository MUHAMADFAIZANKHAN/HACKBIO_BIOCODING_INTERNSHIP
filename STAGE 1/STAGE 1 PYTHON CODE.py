import numpy as np
import pandas as pd

# Task 1: Function for translating DNA to protein
codon_table = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K', 'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*', 'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
}

def translate_dna_to_protein(dna_sequence):
    """
    Translates a DNA sequence into a protein sequence based on the codon table.
    """
    protein_sequence = []
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i+3]
        if codon in codon_table:
            protein_sequence.append(codon_table[codon])
        else:
            protein_sequence.append('X')  # For any codon not in the table
    return ''.join(protein_sequence)

# Example DNA sequence
dna_sequence = "ATGGAGGAGTAA"
protein = translate_dna_to_protein(dna_sequence)
print(f"Protein sequence: {protein}")


# Task 2: Logistic Growth Function

def logistic_growth(t, r, K, L, E):
    """
    Simulates logistic population growth.
    """
    # Ensure L and E are both at least 6 for valid range for randint
    if L <= 5:
        L = 6  # Adjust to be greater than 5
    if E <= 5:
        E = 6  # Adjust to be greater than 5

    # Print values of L and E for debugging
    print(f"L (Lag Phase): {L}, E (Exponential Phase): {E}")
    
    # Random lag and exponential phases
    lag_phase = np.random.randint(5, L)  # Ensures the range is valid
    exp_phase = np.random.randint(5, E)  # Ensures the range is valid

    # Logistic growth curve
    growth = np.zeros(len(t))
    for i in range(len(t)):
        if i < lag_phase:
            growth[i] = 0
        elif i < lag_phase + exp_phase:
            growth[i] = K * (1 - np.exp(-r * (i - lag_phase)))
        else:
            growth[i] = K / (1 + np.exp(-r * (i - lag_phase - exp_phase)))
    
    return growth


# Task 3: Generate a DataFrame with 100 Different Growth Curves

# Time points for simulation
t = np.linspace(0, 100, 100)

# Generate 100 different growth curves with random parameters
growth_data = []
for _ in range(100):
    r = np.random.uniform(0.1, 0.5)  # Growth rate
    K = np.random.uniform(50, 500)   # Carrying capacity
    L = np.random.randint(5, 15)     # Lag phase length
    E = np.random.randint(5, 15)     # Exponential phase length
    growth_curve = logistic_growth(t, r, K, L, E)
    growth_data.append(growth_curve)

# Create DataFrame for the generated curves
df_growth = pd.DataFrame(growth_data).transpose()
df_growth.columns = [f"Curve_{i+1}" for i in range(100)]

# Show a preview of the DataFrame
print(df_growth.head())


# Task 4: Time to Reach 80% of the Maximum Growth (Carrying Capacity)

def time_to_80_percent(K, growth_curve, t):
    """
    Determines the time to reach 80% of the carrying capacity.
    """
    target = 0.8 * K
    # Find the time when the growth curve reaches 80% of carrying capacity
    time_to_reach_target = t[np.where(growth_curve >= target)[0][0]]
    return time_to_reach_target

# Calculate the time to reach 80% of the carrying capacity for the first growth curve
time_to_80 = time_to_80_percent(df_growth.iloc[0, 0], df_growth.iloc[:, 0].values, t)
print(f"Time to reach 80% of carrying capacity: {time_to_80} units")


# Task 5: Hamming Distance between Slack username and Twitter handle

def hamming_distance(str1, str2):
    """
    Calculates the Hamming distance between two strings by comparing their characters at each position.
    """
    max_len = max(len(str1), len(str2))
    str1 = str1.ljust(max_len)
    str2 = str2.ljust(max_len)
    return sum(c1 != c2 for c1, c2 in zip(str1, str2))

# Example Hamming distance between Slack and Twitter
slack_username = "FAIZAN"
twitter_handle = "FAIZAN"
distance = hamming_distance(slack_username, twitter_handle)
print(f"Hamming distance between Slack username and Twitter handle: {distance}")
