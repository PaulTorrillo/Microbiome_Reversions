import pandas as pd
import numpy as np  # For sqrt
import scipy.stats as stats
from scipy.stats import norm



def one_proportion_z_test(success, size, population_proportion):
    if size <= 0:
        return None, None
    sample_proportion = success / size
    SE = np.sqrt(population_proportion * (1 - population_proportion) / size)  # Standard error
    z_stat = (sample_proportion - population_proportion) / SE
    p_value = 2 * norm.sf(np.abs(z_stat))
    return z_stat, p_value


# Load the Excel file
df = pd.read_excel('genes_with_COG.xlsx')

# Initialize dictionaries to store the ratios and sums
results = {}

# COG categories dictionary
cog_categories = {
    'J': 'Translation, ribosomal structure and biogenesis',
    'A': 'RNA processing and modification',
    'K': 'Transcription',
    'L': 'Replication, recombination and repair',
    'B': 'Chromatin structure and dynamics',
    'D': 'Cell cycle control, cell division, chromosome partitioning',
    'Y': 'Nuclear structure',
    'V': 'Defense mechanisms',
    'T': 'Signal transduction mechanisms',
    'M': 'Cell wall/membrane/envelope biogenesis',
    'N': 'Cell motility',
    'Z': 'Cytoskeleton',
    'W': 'Extracellular structures',
    'U': 'Intracellular trafficking, secretion, and vesicular transport',
    'O': 'Posttranslational modification, protein turnover, chaperones',
    'C': 'Energy production and conversion',
    'G': 'Carbohydrate transport and metabolism',
    'E': 'Amino acid transport and metabolism',
    'F': 'Nucleotide transport and metabolism',
    'H': 'Coenzyme transport and metabolism',
    'I': 'Lipid transport and metabolism',
    'P': 'Inorganic ion transport and metabolism',
    'Q': 'Secondary metabolites biosynthesis, transport, and catabolism',
    'R': 'General function prediction only',
    'S': 'Function unknown'
}
# For each valid COG letter
for letter, description in cog_categories.items():
    # Find all rows where the COG Category column contains the description
    filtered_df = df[df['COG Category'].str.contains(description, case=False, na=False)]

    if not filtered_df.empty:
        sum_TCA = filtered_df['TCA Count'].sum()
        sum_TTA = filtered_df['TTA Count'].sum()
        sum_serine = filtered_df['Serine Count'].sum()
        sum_leucine = filtered_df['Leucine Count'].sum()
        count = len(filtered_df)

        if sum_TCA + sum_TTA != 0:
            TTA_TCA_count = sum_TTA + sum_TCA
            ser_leu_count = sum_serine + sum_leucine
        else:
            TTA_TCA_count = 0
            ser_leu_count = 0

        if ser_leu_count > 0:
            results[description] = ((sum_leucine * 0.1299 + sum_serine * 0.1420) / ser_leu_count,TTA_TCA_count, ser_leu_count, count,
                                    15*one_proportion_z_test(TTA_TCA_count, ser_leu_count, (
                                                sum_leucine * 0.1299 + sum_serine * 0.1420) / ser_leu_count)[1])
        else:
            results[description] = (None, None, None, None, None)
    else:
        results[description] = (None, None, None, None, None)

# Output the results
for key, (usage, TTA_TCA_count, ser_leu_count, proteins, test_result) in results.items():
    if TTA_TCA_count is not None and proteins >= 50 and test_result < 0.05:
        print(
            f"COG Category: {key}, Expected Usage: {usage}, TTA/TCA Counts: {TTA_TCA_count}, Serine/Leucine Counts: {ser_leu_count}, Proteins: {proteins}, Test Result: {test_result}")
