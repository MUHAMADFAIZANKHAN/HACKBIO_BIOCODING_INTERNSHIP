#HACKBIO  BIOCODING INTERNSHIP STAGE 2 TASKS

#Import libraries

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import numpy as np

# Task Code 2.4

# Load datasets
sift_path = "https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/R/datasets/sift.tsv"
foldX_path = "https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/R/datasets/foldX.tsv"

# Read datasets with appropriate separators
df_sift = pd.read_csv(sift_path, sep=r"\s+")
df_foldX = pd.read_csv(filepath_or_buffer=foldX_path, sep=r"\s+")

# Create unique identifier for proteins and amino acids
df_sift['specific_Protein_aa'] = df_sift['Protein'].str.cat(df_sift['Amino_Acid'], sep='_')
df_foldX['specific_Protein_aa'] = df_foldX['Protein'].str.cat(df_foldX['Amino_Acid'], sep='_')

# Merge datasets on the unique identifier
merged_df = pd.merge(df_sift[["Protein", "Amino_Acid", 'specific_Protein_aa', 'sift_Score']],
                     df_foldX[['specific_Protein_aa', 'foldX_Score']],
                     on='specific_Protein_aa')

# Filter deleterious mutations (SIFT score < 0.05 and FoldX score > 2)
mutations_deleterious_df = merged_df.loc[(merged_df["sift_Score"] < 0.05) & (merged_df["foldX_Score"] > 2)]

# Generate frequency table of amino acids
frequency_table = {}
for row in mutations_deleterious_df["Amino_Acid"]:
    first_aa = row[0]
    if first_aa not in frequency_table:
        frequency_table[first_aa] = 1
    else:
        frequency_table[first_aa] += 1

# Display frequency table in console
print(frequency_table)

# Lists for amino acids and their frequencies
keys = list(frequency_table.keys())
values = list(frequency_table.values())

# Bar plot for amino acid frequencies
plt.figure(figsize=(10, 6))
plt.bar(keys, values, color='skyblue')
plt.xlabel('Amino Acid')
plt.ylabel('Frequency')
plt.title('Frequency of Amino Acids')
plt.savefig('Bar plot of Frequency of Amino Acids.png')
plt.show()

# Pie chart for amino acid frequencies
plt.figure(figsize=(8, 8))
plt.pie(values, labels=keys, autopct='%1.1f%%', startangle=140)
plt.title('Frequency of Amino Acids')
plt.axis('equal')
plt.savefig('Pie chart of Frequency of Amino Acids.png')
plt.show()


## Amino Acid with Highest Impact:

# Glycine (G) is a small and flexible amino acid that significantly impacts protein structure by allowing tight turns and facilitating folding. It is found 1307 times.

## Structural and Functional Properties of Amino Acids with >100 Occurrences:

# Glycine (G): Small and flexible, aids in protein folding.

# Leucine (L): Hydrophobic, stabilizes protein core.

# Isoleucine (I): Hydrophobic, contributes to structural stability.

# Arginine (R): Basic, involved in binding phosphate groups.

# Proline (P): Causes kinks in protein chains, disrupts regular structures.

# Alanine (A): Small, helps in stabilizing protein structure.

# Tyrosine (Y): Polar, involved in signaling and enzyme activity.

# Valine (V): Hydrophobic, stabilizes protein structure.

# Tryptophan (W): Largest amino acid, involved in protein-protein interactions.    

# Task Code 2.6



# URL of the dataset
url_data = "https://gist.githubusercontent.com/stephenturner/806e31fce55a8b7175af/raw/1a507c4c3f9f1baaa3a69187223ff3d3050628d4/results.txt"

# Load the dataset from the URL
gene_expression_df = pd.read_csv(url_data, sep=" ")

# Calculate the -log10(p-value) for visualization
gene_expression_df['-log10_pvalue'] = -np.log10(gene_expression_df['pvalue'])

## Differential expression analysis 
# Initialize the 'regulation' column with 'No'
gene_expression_df['regulation'] = 'No'

# Label upregulated genes (log2FoldChange > 1 and p-value < 0.01)
gene_expression_df.loc[
    (gene_expression_df['log2FoldChange'] > 1) & (gene_expression_df['pvalue'] < 0.01), 'regulation'
] = 'Up'

# Label downregulated genes (log2FoldChange < -1 and p-value < 0.01)
gene_expression_df.loc[
    (gene_expression_df['log2FoldChange'] < -1) & (gene_expression_df['pvalue'] < 0.01), 'regulation'
] = 'Down'

# Color map for scatter plot
color_map = {'Up': 'red', 'Down': 'blue', 'No': 'grey'}

## Volcano plot for differential expression visualization
plt.figure(figsize=(10, 6))
sns.scatterplot(
    data=gene_expression_df, x='log2FoldChange', y='-log10_pvalue', hue='regulation', palette=color_map
)

# Reference lines for significance thresholds
plt.axhline(y=-np.log10(0.01), color='r', linestyle='--', label='p = 0.01')
plt.axvline(x=1, color='r', linestyle='--', label='log2FC = 1')
plt.axvline(x=-1, color='r', linestyle='--', label='log2FC = -1')

# Labels and title
plt.title('Volcano Plot')
plt.xlabel('log2(Fold Change)')
plt.ylabel('-log10(p-value)')
plt.legend(title='Gene Regulation')
plt.show()

## Extraction of differentially expressed genes
# Upregulated genes
upregulated_genes_df = gene_expression_df.loc[
    (gene_expression_df['log2FoldChange'] > 1) & (gene_expression_df['pvalue'] < 0.01)
]
upregulated_genes_list = upregulated_genes_df["Gene"].tolist()

# Downregulated genes
downregulated_genes_df = gene_expression_df.loc[
    (gene_expression_df['log2FoldChange'] < -1) & (gene_expression_df['pvalue'] < 0.01)
]
downregulated_genes_list = downregulated_genes_df["Gene"].tolist()

## Identification of the most differentially expressed genes
# Sort upregulated genes by log2FoldChange in descending order
upregulated_genes_sorted = upregulated_genes_df.sort_values(by='log2FoldChange', ascending=False)
print("Top 5 upregulated genes:", upregulated_genes_sorted["Gene"].head(5).tolist())

# Sort downregulated genes by log2FoldChange in descending order
downregulated_genes_sorted = downregulated_genes_df.sort_values(by='log2FoldChange', ascending=False)
print("Top 5 downregulated genes:", downregulated_genes_sorted["Gene"].head(5).tolist())

## Top 5 Upregulated Genes:
# 1. DTHD1: This gene encodes a protein containing a death domain, which is involved in signaling pathways and the apoptosis pathway.
# 2. EMILIN2: Encodes a protein that is part of the extracellular matrix and is involved in processes like angiogenesis and platelet aggregation.
# 3. PI16: Encodes a peptidase inhibitor that may inhibit cardiomyocyte growth and is involved in negative regulation of peptidase activity.
# 4. C4orf45: This gene encodes a protein predicted to be localized in the cytoplasm and nucleus, but its specific function is not well characterized.
# 5. FAM180B: Predicted to be located in the extracellular region and is associated with diseases like borderline leprosy and mosaic variegated aneuploidy syndrome.


## Top 5 Downregulated Genes:
# 1. FAM46B: Encodes a polyadenylation polymerase involved in the maintenance of translational efficiency and is essential for the viability of human embryonic stem cells.
# 2. HS3ST3A1: Encodes a sulfotransferase involved in the biosynthesis of heparan sulfate, which plays a role in various biological activities.
# 3. PMEL: Encodes a protein involved in melanosome maturation and melanin production, playing a crucial role in pigmentation.
# 4. TNFAIP6: Encodes a protein involved in extracellular matrix stability and cell migration, and is associated with inflammation.
# 5. COL4A2: Encodes a type IV collagen protein that is a major structural component of basement membranes and is involved in angiogenesis and tumor growth inhibition.

# Task Code 2.7

# Read the CSV file
data_path = "https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/R/nhanes.csv"
df = pd.read_csv(data_path, sep=",")

## Task number 1
# Process all NA
df_process = df.dropna()


## Task number 2
# Create and save the first histogram of BMI
df_BMI = df.dropna(subset=["BMI"])  # Removes rows with missing BMI values from the DataFrame
plt.figure(figsize=(10, 6))  # Sets the figure size for the histogram plot
sns.histplot(df_BMI["BMI"])  # Plots the BMI values using Seaborn's histplot
plt.title('Histogram of BMI')  # Adds a title to the plot
plt.xlabel('BMI')  # Labels the X-axis
plt.savefig('histogram_bmi.png')  # Saves the generated histogram to a PNG file
plt.show()  # Displays the histogram

# Create and save the second histogram of Weight
df_Weight = df.dropna(subset=["Weight"])  # Removes rows with missing Weight values from the DataFrame
plt.figure(figsize=(10, 6))  # Sets the figure size for the histogram plot
sns.histplot(df_Weight["Weight"])  # Plots the Weight values using Seaborn's histplot
plt.title('Histogram of Weight')  # Adds a title to the plot
plt.xlabel('Weight')  # Labels the X-axis
plt.savefig('histogram_weight.png')  # Saves the generated histogram to a PNG file
plt.show()  # Displays the histogram

# Create and save the third histogram of Age
df_Age = df.dropna(subset=["Age"])  # Removes rows with missing Age values from the DataFrame
plt.figure(figsize=(10, 6))  # Sets the figure size for the histogram plot
sns.histplot(df_Age["Age"])  # Plots the Age values using Seaborn's histplot
plt.title('Histogram of Age')  # Adds a title to the plot
plt.xlabel('Age')  # Labels the X-axis
plt.savefig('histogram_age.png')  # Saves the generated histogram to a PNG file
plt.show()  # Displays the histogram

# Create and save the fourth histogram of Weight in pounds
df_Weight_in_Pounds = df.dropna(subset=["Weight"])  # Removes rows with missing Weight values from the DataFrame
df_Weight_in_Pounds = df["Weight"] * 2.2  # Converts Weight from kilograms to pounds
plt.figure(figsize=(10, 6))  # Sets the figure size for the histogram plot
sns.histplot(df_Weight_in_Pounds, bins=20)  # Plots the Weight in pounds values using Seaborn's histplot
plt.title('Histogram of Weight in Pounds')  # Adds a title to the plot
plt.xlabel('Weight in Pounds')  # Labels the X-axis
plt.savefig('histogram_weight_in_pounds.png')  # Saves the generated histogram to a PNG file
plt.show()  # Displays the histogram

## Task number 3

# Calculate the mean of the 'Pulse' column
mean_pulse = df["Pulse"].mean()

# Round the mean value to five decimal places
mean_pulse_rounded = round(mean_pulse, 5)

# Print the rounded mean value
print(f"Mean Pulse (rounded to 5 decimal places): {mean_pulse_rounded}")
# Output: Mean Pulse (rounded to 5 decimal places): 73.63382

## Exercise number 3

# Calculate the maximum value of the 'BPDia' column
max_BPDia = df["BPDia"].max()

# Calculate the minimum value of the 'BPDia' column
min_BPDia = df["BPDia"].min()


# Print the minimum and maximum values
print(f"Minimum Diastolic Blood Pressure: {min_BPDia}")
# Output: Minimum Diastolic Blood Pressure: 0.0
print(f"Maximum Diastolic Blood Pressure: {max_BPDia}")
# Output: Maximum Diastolic Blood Pressure: 116.0

## Task number 4

# Calculate descriptive statistics for the 'Income' column
describe_income = df["Income"].describe()

# Extract the standard deviation (std) from the descriptive statistics
income_std = describe_income['std']

# Calculate the variance directly from the 'Income' column
income_variance = df["Income"].var()

# Round the standard deviation to five decimal places
income_std_rounded = round(income_std, 5)

# Round the variance to five decimal places
income_variance_rounded = round(income_variance, 5)

# Print the rounded standard deviation and variance
print(f"Standard Deviation (std): {income_std_rounded}")
# Output: Standard Deviation (std): 33489.76064
print(f"Variance: {income_variance_rounded}")
print(f"Variance: {income_variance_rounded}")
# Output: Variance: 1121564067.88888

## Task number 5
plt.figure(figsize=(10, 6))  # Set the figure size
sns.scatterplot(data=df, x="Weight", y="Height", hue="Gender")  # Scatter plot colored by Gender
plt.title('Scatterplot of Weight vs Height (Gender)')  # Add a title to the plot
plt.xlabel('Weight')  # Label the X-axis
plt.ylabel('Height')  # Label the Y-axis
plt.savefig('scatterplot_weight_vs_height_gender.png')  # Save the plot as a PNG file
plt.show()  # Display the plot

# Scatterplot 2: Weight vs Height colored by SmokingStatus
plt.figure(figsize=(10, 6))  # Set the figure size
sns.scatterplot(data=df, x="Weight", y="Height", hue="SmokingStatus")  # Scatter plot colored by Smoking Status
plt.title('Scatterplot of Weight vs Height (SmokingStatus)')  # Add a title to the plot
plt.xlabel('Weight')  # Label the X-axis
plt.ylabel('Height')  # Label the Y-axis
plt.savefig('scatterplot_weight_vs_height_smokingstatus.png')  # Save the plot as a PNG file
plt.show()  # Display the plot

# Scatterplot 3: Weight vs Height colored by Diabetes
plt.figure(figsize=(10, 6))  # Set the figure size
sns.scatterplot(data=df, x="Weight", y="Height", hue="Diabetes")  # Scatter plot colored by Diabetes status
plt.title('Scatterplot of Weight vs Height (Diabetes)')  # Add a title to the plot
plt.xlabel('Weight')  # Label the X-axis
plt.ylabel('Height')  # Label the Y-axis
plt.savefig('scatterplot_weight_vs_height_diabetes.png')  # Save the plot as a PNG file
plt.show()  # Display the plot

## Task number 6
# Age vs. Gender
male_age = df[df['Gender'] == 'male']['Age']
female_age = df[df['Gender'] == 'female']['Age']
t_stat, p_value = stats.ttest_ind(male_age.dropna(), female_age.dropna())
print(f"T-Test Age vs. Gender - p-value: {p_value:.5f}")
# Output: T-Test Age vs. Gender - p-value: 0.08020

# BMI vs. Diabetes
bmi_no_diabetes = df[df['Diabetes'] == 'No']['BMI']
bmi_diabetes = df[df['Diabetes'] == 'Yes']['BMI']
t_stat, p_value = stats.ttest_ind(bmi_no_diabetes.dropna(), bmi_diabetes.dropna(), equal_var=False)
print(f"T-Test BMI vs. Diabetes - p-value: {p_value:.5f}")
# Output: T-Test BMI vs. Diabetes - p-value: 0.00000

# Alcohol per year vs. Marital Status
alcohol_single = df[df['RelationshipStatus'] == 'Single']["AlcoholYear"]
alcohol_committed = df[df['RelationshipStatus'] == 'Committed']["AlcoholYear"]
t_stat, p_value = stats.ttest_ind(alcohol_single.dropna(), alcohol_committed.dropna(), equal_var=False)
print(f"T-Test Alcohol vs. Marital Status - p-value: {p_value:.5f}")
# Output: T-Test Alcohol vs. Marital Status - p-value: 0.00000