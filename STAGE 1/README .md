# Hackbio Bio CodingInternship

## Contributors

| Slack Username | Email |
|---------------|-----------------------------|
| @FAIZAN      | faizanjeee@hotmail.com      |
| @Elsaeed     | saiduabdulkadir38@gmail.com |

## Stage 1 Tasks

### Task 1: Function for Translating DNA to Protein

This script contains a function `translate_dna_to_protein()` that translates a DNA sequence into a protein sequence using a predefined codon table. It takes a DNA sequence as input, processes it in triplets (codons), and maps each codon to its corresponding amino acid. If a codon is not found in the table, it is replaced with 'X'.

Example Usage:
```python
# Example DNA sequence
dna_sequence = "ATGGAGGAGTAA"
protein = translate_dna_to_protein(dna_sequence)
print(f"Protein sequence: {protein}")
```

### Task 2: Logistic Growth Function

The script defines a logistic growth function `logistic_growth()` that simulates population growth based on parameters such as growth rate (r), carrying capacity (K), lag phase (L), and exponential phase (E). The function ensures that L and E are valid values for realistic simulations.

### Task 3: Generate a DataFrame with 100 Different Growth Curves

A dataset containing 100 different growth curves is generated using the `logistic_growth()` function. Random values are assigned to the parameters r, K, L, and E to create diverse growth patterns. The resulting data is stored in a Pandas DataFrame and displayed.

### Task 4: Time to Reach 80% of the Maximum Growth (Carrying Capacity)

The function `time_to_80_percent()` calculates the time required for the growth curve to reach 80% of its carrying capacity. It finds the first instance in the dataset where the population reaches this threshold and returns the corresponding time point.

### Task 5: Hamming Distance Calculation

The script includes a `hamming_distance()` function that computes the Hamming distance between two strings by comparing their characters at each position. This is used to compare Slack usernames and Twitter handles.

Example Usage:
```python
slack_username = "FAIZAN"
twitter_handle = "FAIZAN"
distance = hamming_distance(slack_username, twitter_handle)
print(f"Hamming distance between Slack username and Twitter handle: {distance}")
```

This function helps measure the similarity between two strings based on character differences.
