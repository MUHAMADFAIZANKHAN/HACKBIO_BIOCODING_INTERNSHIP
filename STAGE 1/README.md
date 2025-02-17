# Hackbio Bio CodingInternship

## Contributors

| Slack Username | Email |
|---------------|-----------------------------|
| @FAIZAN      | faizanjeee@hotmail.com      |
| @Elsaeed     | saiduabdulkadir38@gmail.com |

## Stage 1 Tasks

### Task 1: Function for Translating DNA to Protein

This script contains a function `translate_dna_to_protein()` that translates a DNA sequence into a protein sequence using a predefined codon table.

#### Code:
```python
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
    protein_sequence = []
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i+3]
        protein_sequence.append(codon_table.get(codon, 'X'))  # 'X' for unknown codons
    return ''.join(protein_sequence)

# Example DNA sequence
dna_sequence = "ATGGAGGAGTAA"
protein = translate_dna_to_protein(dna_sequence)
print(f"Protein sequence: {protein}")
```

### Task 2: Logistic Growth Function

This function simulates logistic population growth based on parameters such as growth rate, carrying capacity, lag phase, and exponential phase.

#### Code:
```python
import numpy as np

def logistic_growth(t, r, K, L, E):
    if L <= 5:
        L = 6
    if E <= 5:
        E = 6

    lag_phase = np.random.randint(5, L)
    exp_phase = np.random.randint(5, E)

    growth = np.zeros(len(t))
    for i in range(len(t)):
        if i < lag_phase:
            growth[i] = 0
        elif i < lag_phase + exp_phase:
            growth[i] = K * (1 - np.exp(-r * (i - lag_phase)))
        else:
            growth[i] = K / (1 + np.exp(-r * (i - lag_phase - exp_phase)))
    
    return growth
```

### Task 3: Generate a DataFrame with 100 Different Growth Curves

A dataset containing 100 different growth curves is generated and stored in a Pandas DataFrame.

#### Code:
```python
import pandas as pd

t = np.linspace(0, 100, 100)
growth_data = []
for _ in range(100):
    r = np.random.uniform(0.1, 0.5)
    K = np.random.uniform(50, 500)
    L = np.random.randint(5, 15)
    E = np.random.randint(5, 15)
    growth_curve = logistic_growth(t, r, K, L, E)
    growth_data.append(growth_curve)

df_growth = pd.DataFrame(growth_data).transpose()
df_growth.columns = [f"Curve_{i+1}" for i in range(100)]
print(df_growth.head())
```

### Task 4: Time to Reach 80% of the Maximum Growth (Carrying Capacity)

This function calculates the time required for a population to reach 80% of its carrying capacity.

#### Code:
```python
def time_to_80_percent(K, growth_curve, t):
    target = 0.8 * K
    time_to_reach_target = t[np.where(growth_curve >= target)[0][0]]
    return time_to_reach_target

time_to_80 = time_to_80_percent(df_growth.iloc[0, 0], df_growth.iloc[:, 0].values, t)
print(f"Time to reach 80% of carrying capacity: {time_to_80} units")
```

### Task 5: Hamming Distance Calculation

Calculates the Hamming distance between two strings by comparing their characters at each position.

#### Code:
```python
def hamming_distance(str1, str2):
    max_len = max(len(str1), len(str2))
    str1 = str1.ljust(max_len)
    str2 = str2.ljust(max_len)
    return sum(c1 != c2 for c1, c2 in zip(str1, str2))

slack_username = "FAIZAN"
twitter_handle = "FAIZAN"
distance = hamming_distance(slack_username, twitter_handle)
print(f"Hamming distance between Slack username and Twitter handle: {distance}")
```
