import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# Configuring plot parameters
plt.rcParams['font.family'] = 'Helvetica'

fig = plt.figure(figsize=(8.5, 8))

# Add first subplot on the top
ax1 = fig.add_subplot(2, 2, 1)

# Add second subplot on the top
ax2 = fig.add_subplot(2, 2, 2)

# Add a larger subplot at the bottom
ax3 = fig.add_subplot(2, 1, 2)

# Defining constants
mutation_rate_per_bp_per_generation = 10**-9
mutation_rate_per_codon_per_generation = 3 * mutation_rate_per_bp_per_generation
s = 3.5 * 10**-5
codons_per_genome = 10**6
core_genome_codons = 500000
fraction_synonymous = 0.25
fraction_nonsynonymous = 0.75
fraction_N_neutral = 1/10

# Calculating rates
neutral_N_rate_per_codon = mutation_rate_per_codon_per_generation * fraction_nonsynonymous * fraction_N_neutral
nonneutral_N_rate_per_codon = mutation_rate_per_codon_per_generation * fraction_nonsynonymous * (1 - fraction_N_neutral)
S_rate_per_codon = mutation_rate_per_codon_per_generation * fraction_synonymous
normalization_constant = fraction_synonymous / fraction_nonsynonymous

# Additional parameters
numpoints = 50000
max_generations = 5000000

# Load data and calculate years
num_points = 10000
x_days = np.arange(num_points) * 10
y_data = np.loadtxt('continuous.txt')
x_years = x_days / 365

# Plot first set of data
ax1.plot(x_days, np.exp(-(0.0010125 * 10000 * (1 - np.exp(-0.35 * 0.0001 * x_days)) / 0.35)),
         linewidth=3, label='theory', color='k', zorder=100)
ax1.plot(x_days, y_data, linewidth=3, label='continuous selection', linestyle='--', color='tab:gray')

# Plot second set of data
ydata2 = np.loadtxt('periodic.txt')
ax1.plot(x_days, ydata2, linewidth=3, linestyle=':', label='periodic selection', color='tab:brown')

# Setting plot attributes
ax1.legend(fontsize=10, frameon=False)
ax1.set(xlabel='bacterial generations', ylabel='frequency of wild type (mutation free) class', yscale='log')
ax1.set_xlabel('bacterial generations', size=14)
ax1.set_ylabel('frequency of wild type (mutation free) class', size=14, rotation='vertical')
ax1.tick_params(length=12, width=1, which='major', direction='inout')
ax1.tick_params(length=8, width=1, which='minor', direction='inout')

# Adding secondary x-axis
ax12 = ax1.twiny()
ax12.set_xlim(ax1.get_xlim())
ax12.set_xlabel('years', size=14, labelpad=10)

# Formatting x-axis labels
def years_formatter1(x, pos):
    return f'{int(x / 365)}'
ax12.xaxis.set_major_formatter(FuncFormatter(years_formatter1))
ax12.tick_params(length=12, width=1, which='major', direction='inout')
ax12.tick_params(length=8, width=1, which='minor', direction='inout')

def small_round(number):
    """Rounds a number up until it reaches a minimum of 100."""
    counter = 0
    rounded_number = number
    while rounded_number < 100:
        rounded_number *= 10
        counter += 1
    rounded_number = np.round(rounded_number)
    rounded_number /= (10 ** counter)
    return rounded_number

def years_formatter(years, pos):
    """Formats years with a scientific notation."""
    formatted_years = small_round((10 ** years) / (2 * S_rate_per_codon * 365))
    return f'{formatted_years:.1e}'

# Generate evenly spaced points between 1 and max_generations
T_generations = np.linspace(1, max_generations, numpoints)

# Initialize empty arrays to hold mutation data
synonymous_muts = np.zeros(numpoints)
neutral_nonsynonymous_muts = np.zeros(numpoints)
dN_dS_theory = np.zeros(numpoints)
theory_fitness = np.zeros(numpoints)

# Load mutation data
mutation_data = [np.loadtxt(f"mutation{i}.txt") for i in range(4)]
dN_dS_data = [np.zeros(numpoints) for _ in range(4)]
fitness_data = [np.ones(numpoints) for _ in range(4)]

# Calculate mutation values and fitness
for n in range(numpoints):
    synonymous_muts[n] = S_rate_per_codon * core_genome_codons * 2 * T_generations[n]
    neutral_nonsynonymous_muts[n] = neutral_N_rate_per_codon * core_genome_codons * 2 * T_generations[n]
    dN_dS_theory[n] = normalization_constant * (
        (neutral_nonsynonymous_muts[n] + 2 * nonneutral_N_rate_per_codon * core_genome_codons * (
            (1 - np.exp(-s * T_generations[n])) / s)) / synonymous_muts[n])
    theory_fitness[n] = (1 - s) ** (
                nonneutral_N_rate_per_codon * core_genome_codons * ((1 - np.exp(-s * T_generations[n])) / s))

    for i in range(4):
        dN_dS_data[i][n] = normalization_constant * (
            (neutral_nonsynonymous_muts[n] + 2 * mutation_data[i][n]) / synonymous_muts[n])
        fitness_data[i][n] = (1 - s) ** mutation_data[i][n]

# Plotting dN/dS ratios for different mutation rates
mutation_labels = ['theory', 'N$_{e}$ 1 million', 'N$_{e}$ 1 billion', 'N$_{e}$ 1 trillion', 'N$_{e}$ 1 quadrillion']
mutations_data = [dN_dS_theory, dN_dS_data[0], dN_dS_data[1], dN_dS_data[2], dN_dS_data[3]]

for label, data in zip(mutation_labels, mutations_data):
    ax3.plot(S_rate_per_codon*2*T_generations, data, linewidth=3, label=label)

# Plot a line at dN/dS = 1 for reference
ax3.plot(S_rate_per_codon*2*T_generations, np.ones_like(dN_dS_theory), linewidth=3, label='d$_N$/d$_S$=1', linestyle='--', color='k')

# Adjust the scales and tick parameters of the axes
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.tick_params(length=12, width=1, which='major', direction='inout')
ax3.tick_params(length=8, width=1, which='minor', direction='inout')

# Set labels for the axes
ax3.set_xlabel('d$_S$ (core genome synonymous divergence)', size=14)
ax3.set_ylabel('d$_N$/d$_S$', size=14, rotation='vertical')

# Set the limits for the axes
ax3.set_xlim([0.00001,0.01])
ax3.set_ylim([0.05,5])

# Display the legend for the plot
ax3.legend(fontsize='10', ncol=2, frameon=False)

# Path to the folder containing flat files
folder_path = "dnds_flat_files"

# Generate a sorted list of unique genera from the file names
genus_list = sorted({filename.split("_")[0] + '_' + filename.split("_")[1] for filename in os.listdir(folder_path) if filename != '.DS_Store'})

# Iterate through the genera list
for genus in genus_list:
    # Load data for each genus into a list
    genus_data = [np.loadtxt(f'dnds_flat_files/{file}', skiprows=1, delimiter=',', usecols=(0, 1))
                  for file in os.listdir(folder_path) if file.startswith(genus)]

    # Concatenate all the data for the current genus
    concatenated_data = np.concatenate(genus_data)

    # Filter data based on threshold and scatter plot on ax3
    if np.sort(concatenated_data[:, 0])[9] < 0.0005:
        x = np.log(concatenated_data[:, 0])
        y = np.log(concatenated_data[:, 1])
        ax3.scatter(np.exp(x), np.exp(y), s=5, alpha=0.2, color='gray')

# Setting the upper X-axis of the plot (ax3)
ax32 = ax3.twiny()
ax32.set_xlabel('MRCA (years)', size=14, labelpad=10)
ax32.set_xlim([-5, -2])
ax32.xaxis.set_major_formatter(FuncFormatter(years_formatter))
ax32.tick_params(length=12, width=1, which='major', direction='inout')
ax32.tick_params(length=8, width=1, which='minor', direction='inout')

# Plotting fitness over time for different effective population sizes
fitnesses = [fitness_data[0], fitness_data[1], fitness_data[2], fitness_data[3]]
labels = ['N$_{e}$ 1 million', 'N$_{e}$ 1 billion', 'N$_{e}$ 1 trillion', 'N$_{e}$ 1 quadrillion']

ax2.plot(T_generations, theory_fitness, linewidth=3, label=label, color='k')

for fitness, label in zip(fitnesses, labels):
    ax2.plot(T_generations, fitness, linewidth=3, label=label)

# Adjust the scales, tick parameters, labels, and legend of ax2
ax2.set_xlabel('bacterial generations', size=14)
ax2.set_ylabel('relative fitness', size=14, rotation='vertical')
ax2.tick_params(length=12, width=1, which='major', direction='inout')
ax2.tick_params(length=8, width=1, which='minor', direction='inout')
ax2.legend(fontsize=10, frameon=False)
ax2.set_xlim([1, 500000])
ax2.set_ylim([0.992, 1])

# Setting the upper X-axis of the plot (ax2)
ax22 = ax2.twiny()
ax22.set_xlim(ax2.get_xlim())
ax22.set_xlabel('years', size=14, labelpad=10)
ax22.xaxis.set_major_formatter(FuncFormatter(years_formatter1))
ax22.tick_params(length=12, width=1, which='major', direction='inout')
ax22.tick_params(length=8, width=1, which='minor', direction='inout')

# Adjusting layout and showing the plot
fig.tight_layout()
plt.show()
