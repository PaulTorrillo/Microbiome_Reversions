import matplotlib
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# Set font size globally for all text elements
matplotlib.rcParams['font.size'] = 10
matplotlib.rcParams['font.family'] = 'Helvetica'


file_path = 'compensatory_mutations.xlsx'
data = pd.read_excel(file_path)

#Use other columns for increased mutation rates
# Pivot the table
pivot_table_corrected = data.pivot(index='Compensatory Rate', columns='Compensatory Weakness', values='Percentage')

# Re-plot the heatmap with the corrected pivot table
plt.figure(figsize=(4, 4))
heatmap_corrected = sns.heatmap(1-pivot_table_corrected, annot=True, fmt=".2f", cmap="YlGnBu",
                                cbar_kws={'label': 'Percentage', 'orientation': 'horizontal'}, vmin=0.45, vmax=1)

# Set title and labels
#plt.title('Heatmap of Output Values by Compensatory Rate and Fitness (Corrected)', fontsize=12)
plt.xlabel('')
plt.ylabel('')

# Adjust the orientation of the labels to make them horizontal
plt.yticks(rotation='horizontal')
# Rotate the color bar label to horizontal
cbar = heatmap_corrected.collections[0].colorbar
cbar.set_label('', rotation=0, horizontalalignment='left')

plt.show()

