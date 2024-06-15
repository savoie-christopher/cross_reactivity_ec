# Import data from CrossDome
import os
import pandas as pd

#folder = "H2DB"
folder = "H2KB"


def get_self_significant(df, threshold):
    return [x for x,y in zip(df["subject"], df["pvalue"]) if y < threshold]


# Dictionary of crossdome files.
crossdome_peptides = {}
for peptide_file in os.listdir(f"{folder}/"):
    crossdome_peptides[peptide_file.replace(".csv", "")] = pd.read_csv(f"{folder}/{peptide_file}")
    
# Get list of all cross-acting self-peptides
self_peptides = []
for peptide, df in crossdome_peptides.items():
    self_peptides.append(x for x in get_self_significant(df, 0.0001))

self_peptides = [item for sublist in self_peptides for item in sublist]
self_peptides = list(set(self_peptides))


# Define a function to calculate the value based on indices
def calculate_value(i, j):
    df = crossdome_peptides[i]
    return -np.log10(list(df[df["subject"]==j]["pvalue"])[0])
    
    
x = len(crossdome_peptides)
y = len(self_peptides)

# Create a DataFrame using a dictionary comprehension
data = {f'{j}': [calculate_value(i, j) for i in crossdome_peptides.keys()] for j in self_peptides}
df = pd.DataFrame(data, index=[f'{i}' for i in crossdome_peptides.keys()])
df.sort_index(inplace=True)


#%%
# Make plot heatmap of p-values
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap

max_value = df.to_numpy().max()

# Define a function to create a custom colormap
def create_custom_colormap(threshold, low_color, high_color):
    colors = [(low_color), (low_color), (high_color)]
    nodes = [0.0, threshold / max_value, 1.0]
    cmap = LinearSegmentedColormap.from_list("custom_cmap", list(zip(nodes, colors)))
    return cmap

# Define the threshold and colors
threshold = 3.99
low_color = "white"
high_color = "red"

# Create the custom colormap
custom_cmap = create_custom_colormap(threshold, low_color, high_color)

# Create the heatmap
fig, ax = plt.subplots(figsize=(10, 8))
cax = ax.imshow(df, cmap=custom_cmap, aspect='auto')
cbar = fig.colorbar(cax, ax=ax, orientation='vertical')
cbar.set_label('-Log10 (p-value)')

# Add title and labels
plt.title('Cross-reactivity of non-self peptides')
plt.xlabel("Self-peptides presented across different cell types and tissues (H2KB-9mers)")
plt.ylabel("Non-self peptides (H2KB-9mers)")

plt.xticks(ticks=np.arange(y), labels=list(df.columns), rotation=90)
plt.yticks(ticks=np.arange(x), labels=df.index)

# Save the figure with higher resolution
plt.savefig('Cross-reactivity of non-self peptides-H2KB.png', dpi=900, bbox_inches='tight')

# Show the figure
plt.show()