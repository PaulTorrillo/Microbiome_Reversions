import numpy as np
import matplotlib.pyplot as plt
import os
import scipy
from matplotlib.ticker import FuncFormatter
import matplotlib

# Configuring matplotlib settings
matplotlib.rcParams['font.family'] = 'Helvetica'

# Defining constants
mutation_rate_per_base_pair_per_generation = 10**-9
mutation_rate_per_codon_per_generation = 3 * mutation_rate_per_base_pair_per_generation
codons_per_genome = 10**6
core_genome_codons = 500000
fraction_synonymous = 0.25
fraction_nonsynonymous = 0.75
fraction_N_neutral = 1 / 10
synonymous_rate_per_codon = mutation_rate_per_codon_per_generation * fraction_synonymous

def to_year(x):
    years = x / (2 * synonymous_rate_per_codon * 365)
    return years

def to_dS(x):
    dS=x*(2 * synonymous_rate_per_codon * 365)
    return dS

# Function for calculating log result
def calculate_log_result(x, s, a):
    t = x / (2 * synonymous_rate_per_codon)
    neutral_non_syn = a * mutation_rate_per_codon_per_generation * fraction_nonsynonymous * core_genome_codons * 2 * t
    synonymous = 2 * synonymous_rate_per_codon * core_genome_codons * t
    result = (1 / 3) * ((neutral_non_syn + 2 * (1 - a) * mutation_rate_per_codon_per_generation * fraction_nonsynonymous * core_genome_codons * ((1 - np.exp(-s * t)) / (s))) / (synonymous))
    return np.log(result)

# Function for rounding
def small_round(x):
    counter = 0
    newx = x
    if newx < 100:
        while newx < 10:
            newx = newx * 10
            counter = counter + 1
        newx = np.round(newx)
        newx = newx / (10 ** counter)
    else:
        while newx > 100:
            newx = newx / 10
            counter = counter + 1
        newx = np.round(newx)
        newx = newx * (10 ** counter)
    return newx

# Function for formatting years
def years_formatter(x, pos):
    years = small_round((10 ** x) / (2 * synonymous_rate_per_codon * 365))
    return years

# Function for calculating coefficient
def calculate_coefficient(x, y):
    popt, pcov = scipy.optimize.curve_fit(calculate_log_result, np.exp(x), y, p0=[0.00001, 0.1])
    return popt[0]

# Function for italicizing text
def italicize(text):
    return r'$\it{%s}$' % text

# Specifying folder path
folder_path = "../dnds_flat_files"
genus_list = set()
# Set the colormap to 'tab10'
cmap = plt.get_cmap('tab10')
# Listing all genus
for file in os.listdir(folder_path):
    genus_list.add(file.split("_")[0] + '_' + file.split("_")[1])
genus_list = [*genus_list]
genus_list.sort()



# Adding data to the plot

genus_counter = -1

# Looping over all genus in the list
for genus in genus_list:
    genus_to_combine = []

    # Collect all files corresponding to a genus
    for file in os.listdir(folder_path):
        if file.startswith(genus):
            genus_to_combine.append(np.loadtxt('../dnds_flat_files/' + file, skiprows=1, delimiter=',', usecols=(0, 1)))

    # Combine all files of a genus
    combined_genus = genus_to_combine[0]
    for i in range(len(genus_to_combine)):
        if i > 0:
            combined_genus = np.concatenate((combined_genus, genus_to_combine[i]), axis=0)

    # Calculate and plot if condition met
    if np.sort(combined_genus[:, 0])[9] < 0.0005:
        x = np.log(combined_genus[:, 0])
        y = np.log(combined_genus[:, 1])
        myline = np.linspace(np.log(0.000001), np.log(1), 1000)
        popt, pcov = scipy.optimize.curve_fit(calculate_log_result, np.exp(x), y, p0=[1E-5, 0.1])


        genus_counter = genus_counter + 1
        color = cmap(genus_counter)
        plt.plot(np.exp(myline), np.exp(calculate_log_result(np.exp(myline), *popt)),
                    label=italicize(genus).replace('Phocaeicola_', 'Ph. ').replace('Bacteroides_', 'B. ').replace('Parabacteroides_', 'Pa. '),
                    linewidth='3',color=color)
        plt.scatter(np.exp(x), np.exp(y), s=10, alpha=0.25,color=color)
        plt.xscale('log')
        plt.yscale('log')
        # Set labels for the x and y axes of the main chart
        plt.xlabel('dS (core genome synonymous divergence)', size=14)
        plt.ylabel('dN/dS', size=14, rotation='vertical')

        # Set boundaries for the x and y axes of the main chart
        plt.xlim([0.0000008, 0.1])
        plt.ylim([0.05, 10])
        plt.tick_params(length=10, width=1, which='major', direction='inout', labelsize=12)
        plt.tick_params(length=6, width=1, which='minor', direction='inout')
        plt.legend(fontsize=16)
        plt.tight_layout()
        plt.show()


plt.show()
