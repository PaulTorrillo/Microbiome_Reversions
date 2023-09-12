import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# Setting the font to be used
matplotlib.rcParams['font.family'] = 'Helvetica'

# Defining variables related to mutation rates
mutation_rate_per_bp_per_generation = 10**-9
mutation_rate_per_codon_per_generation = 3 * mutation_rate_per_bp_per_generation
core_genome_codons = 500000
fraction_synonymous = 0.25
fraction_nonsynonymous = 0.75
reversion_penalty=0.2 #Reversions are going to be a specific amino acid mutation so we add a 1 in 5 penalty given mutation is nonsynonymous

# Calculating mutation rates
nonneutral_N_rate_per_codon = mutation_rate_per_codon_per_generation * fraction_nonsynonymous * reversion_penalty

# Initializing arrays for population sizes and calculated properties
population_sizes_logarithmic = np.linspace(3, 15, 1000)  # Population sizes in log scale
prob_fix_revert = np.zeros(1000)  # Probability of fixing revertant per generation
expected_years_to_reversion = np.zeros(1000)  # Expected years to reversion

# Arrays for different selection coefficient scenarios
selection_coefficients = [0.3, 0.03, 0.003]
prob_fix_revert_for_s = {s: np.zeros(1000) for s in selection_coefficients}
expected_years_to_reversion_for_s = {s: np.zeros(1000) for s in selection_coefficients}

# Calculate properties for each population size
for i, population_size_logarithmic in enumerate(population_sizes_logarithmic):
    population_size = 10 ** population_size_logarithmic
    effective_mutation_rate = min(population_size * nonneutral_N_rate_per_codon, 1)
    max_value = max(1, population_size * nonneutral_N_rate_per_codon)

    # Calculate properties for different selection coefficient scenarios
    for s in selection_coefficients:
        prob_fix_revert_for_s[s][i] = effective_mutation_rate * (1 - (1-2*s) ** max_value)
        expected_years_to_reversion_for_s[s][i] = (
            (effective_mutation_rate * (1 - (1-2*s) ** max_value)) ** -1 +
            (np.log(population_size / max_value) / s) / 365
        )

# Create a figure and a set of subplots (2 rows, 1 column)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4, 3.8))
fig.subplots_adjust(left=0.2, right=0.8, bottom=0.2)

# Plotting the probabilities and expected years for different selection coefficient scenarios
line_styles = [':', '--', '-']
lines = []
for s, line_style in zip(selection_coefficients, line_styles):
    line, = ax1.loglog(10 ** population_sizes_logarithmic, prob_fix_revert_for_s[s], label=f's={s}', linewidth='3', color='k', linestyle=line_style)
    lines.append(line)
    ax2.loglog(10 ** population_sizes_logarithmic, expected_years_to_reversion_for_s[s], linewidth='3', color='k', linestyle=line_style)

# Adjust the scales, labels, and tick parameters of the axes
ax1.set_xlim([10 ** 6, 10 ** 12])
ax1.set_ylim([10 ** -6, 10 ** 0])
ax1.set_ylabel('probability fixing revertant \narises per generation', fontsize=10)
ax1.tick_params(axis='both', which='major', labelsize=8)

ax2.set_xlim([10 ** 6, 10 ** 12])
ax2.set_ylim([10 ** 0, 10 ** 3])
ax2.set_xlabel('population size', fontsize=10)
ax2.set_ylabel('expected years to \nfixation of reversion', fontsize=10)
ax2.tick_params(axis='both', which='major', labelsize=8)

# Adjust the layout and add a legend to the figure
fig.subplots_adjust(right=0.7)
labels = [line.get_label() for line in lines]
fig.legend(lines, labels, loc = 'center right', title='selective\nadvantage', fontsize=10, frameon=False)

plt.show()
