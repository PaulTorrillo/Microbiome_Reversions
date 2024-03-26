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
def calculate_result(x, s, a):
    t = x / (2 * synonymous_rate_per_codon)
    neutral_non_syn = a * mutation_rate_per_codon_per_generation * fraction_nonsynonymous * core_genome_codons * 2 * t
    synonymous = 2 * synonymous_rate_per_codon * core_genome_codons * t
    result = (1 / 3) * ((neutral_non_syn + 2 * (1 - a) * mutation_rate_per_codon_per_generation * fraction_nonsynonymous * core_genome_codons * ((1 - np.exp(-s * t)) / (s))) / (synonymous))
    return result

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
    return f'{years:.1e}'

# Function for calculating coefficient
def calculate_coefficient(x, y):
    popt, pcov = scipy.optimize.curve_fit(calculate_result, x, y, p0=[1E-6, 0.1], bounds=([0, 0], [np.inf, np.inf]))
    return popt[0]

# Function for italicizing text
def italicize(text):
    return r'$\it{%s}$' % text

# Specifying folder path
folder_path = "../dnds_flat_files"
genus_list = set()

# Listing all genus
for file in os.listdir(folder_path):
    genus_list.add(file.split("_")[0] + '_' + file.split("_")[1])
genus_list = [*genus_list]
genus_list.sort()

# Plot settings
fig, axs = plt.subplots(2, 1, figsize=(8.5, 8.2), gridspec_kw={'height_ratios': [1.9, 1.1]})
ax_inset = axs[0].inset_axes([0.66, 0.54, 0.3, 0.4])

genus_counter = 1

collected_rs=[]

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
        x = combined_genus[:, 0]
        y = combined_genus[:, 1]
        myline = np.linspace(np.log(0.000001), np.log(1), 1000)
        popt, pcov = scipy.optimize.curve_fit(calculate_result, x, y, p0=[1E-6, 0.1], bounds=([0, 0], [np.inf, np.inf]))
        result = scipy.stats.bootstrap((x, y), calculate_coefficient, paired=True, n_resamples=999, method='basic')
        lower, upper = result.confidence_interval
        lower_error = popt[0] - lower
        lower_error = min(popt[0], lower_error)
        upper_error = upper - popt[0]
        residuals = y - calculate_result(x, *popt)
        var = y - np.mean(y)
        rss = np.sum(residuals ** 2)
        tss = np.sum(var ** 2)
        z = round(1 - (rss / tss), 3)
        collected_rs.append(z)
        axs[1].errorbar(genus_counter, popt[0], yerr=[[lower_error], [upper_error]], fmt='o', capsize=3)
        genus_counter = genus_counter + 1
        axs[0].scatter(x, y, s=10, alpha=0.25)
        ax_inset.scatter(x, y * x, s=10, alpha=0.25)
        axs[0].plot(np.exp(myline), calculate_result(np.exp(myline), *popt),
                    label=italicize(genus).replace('Phocaeicola_', 'Ph. ').replace('Bacteroides_', 'B. ').replace(
                        'Parabacteroides_', 'Pa. '),
                    linewidth='3')
        ax_inset.plot(np.exp(myline), calculate_result(np.exp(myline), *popt) * np.exp(myline))

# Further plot configuration and display

# Set boundaries for the inset chart
ax_inset.set_xlim([0.00001, 0.0003])

# Set axes of main plot to logarithmic scale
axs[0].set_xscale('log')
axs[0].set_yscale('log')

# Set y label for the lower chart
axs[1].set_ylabel('s (selective coefficent)', rotation='vertical', size=14)

# Set boundaries for the y axis of the inset chart
ax_inset.set_ylim([0.00001, 0.0001])

# Set labels for the inset chart axes
ax_inset.set_xlabel('d$_S$', size=12)
ax_inset.set_ylabel('d$_N$', size=12)

# Set tick parameters for the inset and main charts
ax_inset.tick_params(length=8, width=1, which='major', direction='inout', labelsize=10)
axs[0].tick_params(length=10, width=1, which='major', direction='inout', labelsize=12)
axs[0].tick_params(length=6, width=1, which='minor', direction='inout')

# Hide the x ticks for the lower chart and set tick parameters
axs[1].set_xticks([])
axs[1].tick_params(length=10, width=1, which='major', direction='inout', labelsize=12)
axs[1].tick_params(axis='x', which='both', bottom=False, top=False)
axs[1].set_ylim([-1E-5,1E-4])
# Set labels for the x and y axes of the main chart
axs[0].set_xlabel('d$_S$ (core genome synonymous divergence)', size=14)
axs[0].set_ylabel('d$_N$/d$_S$', size=14, rotation='vertical')

# Set boundaries for the x and y axes of the main chart
axs[0].set_xlim([0.0000008, 0.1])
axs[0].set_ylim([0.05, 10])

# Configure legends
lines, labels = axs[0].get_legend_handles_labels()
axs[0].legend(lines[0:],labels[0:], loc='upper center', bbox_to_anchor=(0.5,-0.2),ncol=5,labelspacing=0.5,fontsize=10,frameon=False)

# Configure the secondary x axis of the main chart
ax2 = axs[0].secondary_xaxis('top', functions=(to_year,to_dS))
ax2.set_xlabel('MRCA (years)', size=14)
ax2.tick_params(length=10, width=1, which='major', direction='inout',labelsize=12)
ax2.tick_params(length=6, width=1, which='minor', direction='inout')

print(np.median(collected_rs))
# Ensure a tight layout for the plot and show it
plt.tight_layout()
plt.show()
