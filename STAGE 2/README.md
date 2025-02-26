# HackBio BioCoding Internship - Stage 2 Submission

## Task 2.4: Protein Mutation Analysis

### Objective

This project aims to combine and analyze data from the SIFT and FoldX datasets to investigate protein mutations. The SIFT dataset provides functional impact scores for amino acid substitutions, where a score below 0.05 is considered deleterious. The FoldX dataset provides structural impact scores, with scores above 2 kCal/mol considered deleterious.

We identify mutations that are both functionally and structurally impactful, analyze the frequency of amino acids involved in these mutations, and visualize the data using bar plots and pie charts.

### Steps to Complete the Task

#### 1. Import Datasets

```python
import pandas as pd
import matplotlib.pyplot as plt

# Load datasets
sift_path = "https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/R/datasets/sift.tsv"
foldX_path = "https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/R/datasets/foldX.tsv"
df_sift = pd.read_csv(sift_path, sep=r"\s+")
df_foldX = pd.read_csv(foldX_path, sep=r"\s+")
```

#### 2. Create Unique Identifier

```python
df_sift['specific_Protein_aa'] = df_sift['Protein'].str.cat(df_sift['Amino_Acid'], sep='_')
df_foldX['specific_Protein_aa'] = df_foldX['Protein'].str.cat(df_foldX['Amino_Acid'], sep='_')
```

#### 3. Merge Datasets

```python
merged_df = pd.merge(df_sift[["Protein", "Amino_Acid", 'specific_Protein_aa', 'sift_Score']],
                     df_foldX[['specific_Protein_aa', 'foldX_Score']],
                     on='specific_Protein_aa')
```

#### 4. Identify Deleterious Mutations

```python
mutations_deleterious_df = merged_df.loc[(merged_df["sift_Score"] < 0.05) & (merged_df["foldX_Score"] > 2)]
```

#### 5. Analyze Amino Acid Substitution

```python
frequency_table = {}
for row in mutations_deleterious_df["Amino_Acid"]:
    first_aa = row[0]
    frequency_table[first_aa] = frequency_table.get(first_aa, 0) + 1
```

#### 6. Generate Frequency Table

```python
keys = list(frequency_table.keys())
values = list(frequency_table.values())
```

#### 7. Visualize Data

**Bar Plot:**

```python
plt.figure(figsize=(10, 6))
plt.bar(keys, values, color='skyblue')
plt.xlabel('Amino Acid')
plt.ylabel('Frequency')
plt.title('Frequency of Amino Acids')
plt.savefig('Bar_plot_of_Frequency_of_Amino_Acids.png')
plt.show()
```
(FIGURES/Bar%20plot%20of%20Frequency%20of%20Amino%20Acids.png)

**Pie Chart:**

```python
plt.figure(figsize=(8, 8))
plt.pie(values, labels=keys, autopct='%1.1f%%', startangle=140)
plt.title('Frequency of Amino Acids')
plt.axis('equal')
plt.savefig('Pie_chart_of_Frequency_of_Amino_Acids.png')
plt.show()
```

---

## Task 2.6: Gene Expression Analysis

### Objective

This project analyzes a gene expression dataset, performing data cleaning, visualization, and differential expression analysis. The dataset contains gene expression values, p-values, and log2 fold changes, helping identify upregulated and downregulated genes.

### Steps to Complete the Task

#### 1. Load the Dataset

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load the dataset
url_data = "https://gist.githubusercontent.com/stephenturner/806e31fce55a8b7175af/raw/1a507c4c3f9f1baaa3a69187223ff3d3050628d4/results.txt"
gene_expression_df = pd.read_csv(url_data, sep=" ")
```

#### 2. Compute -log10(p-value)

```python
gene_expression_df['-log10_pvalue'] = -np.log10(gene_expression_df['pvalue'])
```

#### 3. Differential Expression Analysis

```python
gene_expression_df['regulation'] = 'No'
gene_expression_df.loc[(gene_expression_df['log2FoldChange'] > 1) & (gene_expression_df['pvalue'] < 0.01), 'regulation'] = 'Up'
gene_expression_df.loc[(gene_expression_df['log2FoldChange'] < -1) & (gene_expression_df['pvalue'] < 0.01), 'regulation'] = 'Down'
```

#### 4. Volcano Plot Visualization

```python
plt.figure(figsize=(10, 6))
sns.scatterplot(data=gene_expression_df, x='log2FoldChange', y='-log10_pvalue', hue='regulation', palette={'Up': 'red', 'Down': 'blue', 'No': 'grey'})
plt.axhline(y=-np.log10(0.01), color='r', linestyle='--', label='p = 0.01')
plt.axvline(x=1, color='r', linestyle='--', label='log2FC = 1')
plt.axvline(x=-1, color='r', linestyle='--', label='log2FC = -1')
plt.title('Volcano Plot')
plt.xlabel('log2(Fold Change)')
plt.ylabel('-log10(p-value)')
plt.legend(title='Gene Regulation')
plt.show()
```

#### 5. Extracting Differentially Expressed Genes

```python
upregulated_genes_df = gene_expression_df.loc[(gene_expression_df['log2FoldChange'] > 1) & (gene_expression_df['pvalue'] < 0.01)]
upregulated_genes_list = upregulated_genes_df["Gene"].tolist()
downregulated_genes_df = gene_expression_df.loc[(gene_expression_df['log2FoldChange'] < -1) & (gene_expression_df['pvalue'] < 0.01)]
downregulated_genes_list = downregulated_genes_df["Gene"].tolist()
```

---

## Task 2.7: NHANES Dataset Analysis

### Objective

This project analyzes the NHANES dataset, performing data cleaning, visualization, and statistical analysis. The dataset includes health-related attributes such as BMI, weight, age, and lifestyle factors.

### Steps to Complete the Task

#### 1. Load the Dataset

```python
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

data_path = "https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/R/nhanes.csv"
df = pd.read_csv(data_path, sep=",")
```

#### 2. Data Preprocessing

```python
df_process = df.dropna()
```

#### 3. Generate Histograms

```python
plt.figure(figsize=(10, 6))
sns.histplot(df_process["BMI"], kde=True)
plt.title('Histogram of BMI')
plt.xlabel('BMI')
plt.show()
```

---

## Conclusion

The above analyses provide insights into protein mutations, gene expression patterns, and health attributes from the NHANES dataset. These methods allow for data-driven conclusions and visualization of biological trends.
