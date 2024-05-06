import pandas as pd
import numpy as np  # For sqrt
import scipy.stats as stats
from scipy.stats import norm

def one_proportion_z_test(success, size, population_proportion):
    """
    Perform a one-proportion z-test.

    Parameters:
    - success: number of successes in the sample
    - size: sample size
    - population_proportion: the known or hypothesized proportion in the population

    Returns:
    - z_stat: the calculated z-statistic
    - p_value: the p-value associated with the z-statistic
    """
    if size <= 0:
        return None, None

    sample_proportion = success / size
    SE = np.sqrt(population_proportion * (1 - population_proportion) / size)  # Standard error

    # Calculate the z statistic
    z_stat = (sample_proportion - population_proportion) / SE

    # Calculate the p-value
    p_value = 2 * norm.sf(np.abs(z_stat))

    return z_stat, p_value

# Load the Excel file
df = pd.read_excel('genes_with_COG_and_location.xlsx')

# Initialize dictionaries to store the ratios, standard deviations, and sums
results = {}

categories = ['Cytoplasmic', 'CytoplasmicMembrane', 'Extracellular', 'OuterMembrane', 'Periplasmic']

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
    for category in categories:

        # Find all rows where the first column contains the letter and the fourth column matches the category
        filtered_df = df[df['COG Category'].str.contains(description, case=False, na=False) & (df['Cellular Location'] == category)]

        # Check if there are rows to perform the operation
        if not filtered_df.empty:
            # Calculate the sum of all values in column 2 and column 3
            sum_TCA = filtered_df['TTA Count'].sum()
            sum_TTA = filtered_df['TCA Count'].sum()
            sum_serine=filtered_df['Serine Count'].sum()
            sum_lecuine=filtered_df['Leucine Count'].sum()
            count = len(filtered_df['TTA Count'])
            # Avoid division by zero for ratio calculation
            if sum_TCA+sum_TTA != 0:
                TTA_TCA_count = sum_TTA+sum_TCA
                # Calculate the standard deviation
                ser_leu_count = sum_serine+sum_lecuine
            else:
                TTA_TCA_count = 0
                ser_leu_count = 0

            # Include sum_col3 in the results
            combination_key = f"{description}_{category}"
            if ser_leu_count>0:
                #If with both categories
                results[combination_key] = ((sum_lecuine*0.1268+sum_serine*0.1409)/ser_leu_count,ser_leu_count, count, 15*one_proportion_z_test(TTA_TCA_count, ser_leu_count, (sum_lecuine*0.1268+sum_serine*0.1409)/ser_leu_count)[1])
            else:
                combination_key = f"{description}_{category}"
                results[combination_key] = (None, None, None, None)
        else:
            combination_key = f"{description}_{category}"
            results[combination_key] = (None, None, None, None)

# Output the results
for key, (TTA_TCA_count, ser_leu_count,proteins, test_result) in results.items():
    # Check if the ratio is not None before printing
    if TTA_TCA_count is not None and proteins>=50 and test_result<0.05:
        print(f"Combination: {key}, TTA/TCA Counts: {TTA_TCA_count}, Serine/Leucine Counts: {ser_leu_count}, Proteins: {proteins}, Test Result: {test_result}")
