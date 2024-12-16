#bubble chart of figure 3E

import matplotlib.pyplot as plt
from collections import Counter
import pandas as pd
from pandas.api.types import CategoricalDtype

# to read the data file
file_path = '/Users/avanini/riboseq/Fusion_Transcripts_THELAST/data_for_bubble_chart_all.txt'
data = pd.read_csv(file_path, delimiter='\t', names=['Gene1', 'Characteristic1', 'Gene2', 'Characteristic2'])

# to extract and count the pairs of categories
pair_counts = Counter((row['Characteristic1'], row['Characteristic2']) for idx, row in data.iterrows())

# Convert to a DataFrame for easier manipulation
df = pd.DataFrame(pair_counts.items(), columns=["Pair", "Count"])
df["Transcript1"] = df["Pair"].apply(lambda x: x[0])
df["Transcript2"] = df["Pair"].apply(lambda x: x[1])

# Define the desired order of categories
categories_order = [
    'CDS', 'UTR', 'exonic(no-known-CDS)', 'intronic'
]

# Set the order for the categorical data
cat_type = CategoricalDtype(categories=categories_order, ordered=True)
df["Transcript1"] = df["Transcript1"].astype(cat_type)
df["Transcript2"] = df["Transcript2"].astype(cat_type)

# Sort the DataFrame to ensure the order is applied in the plot
df = df.sort_values(by=["Transcript1", "Transcript2"])

# Create the bubble chart
plt.figure(figsize=(10, 10))
plt.scatter(df["Transcript2"], df["Transcript1"], s=df["Count"]*50, color='lightgrey', edgecolors="w", linewidth=0.5)  

# Add labels to bubbles
for i, row in df.iterrows():
    plt.text(row["Transcript2"], row["Transcript1"], row["Count"], ha='center', va='center', fontsize=10)  

# Add labels and title
plt.xlabel('Downstream fusion gene')
plt.ylabel('Upstream fusion gene')
plt.xticks(ticks=range(len(categories_order)), labels=categories_order, rotation=90)
plt.yticks(ticks=range(len(categories_order)), labels=categories_order[::-1])
plt.grid(False)
plt.tight_layout()

# to extend the limits of the axes
plt.xlim(-0.5, len(categories_order) - 0.5)
plt.ylim(len(categories_order) - 0.5, -0.5)

plt.show()
