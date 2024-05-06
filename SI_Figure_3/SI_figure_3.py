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
    return f'{years:.1e}'

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

# Listing all genus
for file in os.listdir(folder_path):
    genus_list.add(file.split("_")[0] + '_' + file.split("_")[1])
genus_list = [*genus_list]
genus_list.sort()

# Plot settings
fig, axs = plt.subplots(2, 1, figsize=(8.5, 8.2), gridspec_kw={'height_ratios': [1, 1]})

# Adding data to the plot


genus_counter = 1

x_mega=[]
y_mega=[]
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
        x_mega=np.append(x_mega,x)
        y = np.log(combined_genus[:, 1])
        y_mega = np.append(y_mega, y)
        myline = np.linspace(np.log(0.000001), np.log(1), 1000)
        popt, pcov = scipy.optimize.curve_fit(calculate_log_result, np.exp(x), y, p0=[1E-5, 0.1])
        genus_counter = genus_counter + 1
        axs[0].scatter(np.exp(x), np.exp(y), s=10, alpha=0.25,color='gray')
        axs[1].scatter(np.exp(x), np.exp(y) * np.exp(x), s=10, alpha=0.25,color='gray')
        #axs[0].plot(np.exp(myline), np.exp(calculate_log_result(np.exp(myline), *popt)),
            #        label=italicize(genus).replace('Bacteroides_', 'B. ').replace('Parabacteroides_', 'P. '),
           #         linewidth='3')
        #axs[1].plot(np.exp(myline), np.exp(calculate_log_result(np.exp(myline), *popt)) * np.exp(myline))

# Further plot configuration and display
residuals=y_mega-calculate_log_result(np.exp(x_mega),3.5E-5,0.10)
var=y_mega-np.mean(y_mega)
rss=np.sum(residuals**2)
tss = np.sum(var ** 2)
z=round(1-(rss/tss),3)

axs[0].plot(np.exp(myline), np.exp(calculate_log_result(np.exp(myline), 3.5E-5,0.1)),
               label='s=3.5 x 10$^{-5}$, R$^2$ = '+str(z),
              linewidth='3')
axs[1].plot(np.exp(myline), np.exp(calculate_log_result(np.exp(myline), 3.5E-5,0.1)) * np.exp(myline),
               label='s=3.5 x 10$^{-5}$, R$^2$ = '+str(z),
              linewidth='3')
residuals=y_mega-calculate_log_result(np.exp(x_mega),3.5E-4,0.1)
var=y_mega-np.mean(y_mega)
rss=np.sum(residuals**2)
tss = np.sum(var ** 2)
z=round(1-(rss/tss),3)
print(z)
axs[0].plot(np.exp(myline), np.exp(calculate_log_result(np.exp(myline), 3.5E-4,0.1)),
               label='s=3.5 x 10$^{-4}$, R$^2$ = '+str(z),
              linewidth='3')
axs[1].plot(np.exp(myline), np.exp(calculate_log_result(np.exp(myline), 3.5E-4,0.1)) * np.exp(myline),
               label='s=3.5 x 10$^{-4}$, R$^2$ = '+str(z),
              linewidth='3')
residuals=y_mega-calculate_log_result(np.exp(x_mega),3.5E-3,0.1)
var=y_mega-np.mean(y_mega)
rss=np.sum(residuals**2)
tss = np.sum(var ** 2)
z=round(1-(rss/tss),3)

axs[0].plot(np.exp(myline), np.exp(calculate_log_result(np.exp(myline), 3.5E-3,0.1)),
               label='s=3.5 x 10$^{-3}$, R$^2$ = '+str(z),
              linewidth='3')
axs[1].plot(np.exp(myline), np.exp(calculate_log_result(np.exp(myline), 3.5E-3,0.1)) * np.exp(myline),
               label='s=3.5 x 10$^{-3}$, R$^2$ = '+str(z),
              linewidth='3')
residuals=y_mega-calculate_log_result(np.exp(x_mega),3.5E-2,0.1)
var=y_mega-np.mean(y_mega)
rss=np.sum(residuals**2)
tss = np.sum(var ** 2)
z=round(1-(rss/tss),3)

axs[0].plot(np.exp(myline), np.exp(calculate_log_result(np.exp(myline), 3.5E-2,0.1)),
               label='s=3.5 x 10$^{-2}$, R$^2$ = '+str(z),
              linewidth='3')
axs[1].plot(np.exp(myline), np.exp(calculate_log_result(np.exp(myline), 3.5E-2,0.1)) * np.exp(myline),
               label='s=3.5 x 10$^{-2}$, R$^2$ = '+str(z),
              linewidth='3')
# Set boundaries for the inset chart
axs[1].set_xlim([0.0000, 0.0006])

# Set axes of main plot to logarithmic scale
axs[0].set_xscale('log')
axs[0].set_yscale('log')



# Set boundaries for the y axis of the inset chart
axs[1].set_ylim([0.0000, 0.0002])

# Set labels for the inset chart axes
axs[1].set_xlabel('d$_S$', size=12)
axs[1].set_ylabel('d$_N$', size=12)

# Set tick parameters for the inset and main charts
axs[1].tick_params(length=8, width=1, which='major', direction='inout', labelsize=10)
axs[0].tick_params(length=10, width=1, which='major', direction='inout', labelsize=12)
axs[0].tick_params(length=6, width=1, which='minor', direction='inout')


# Set labels for the x and y axes of the main chart
axs[0].set_xlabel('d$_S$ (core genome synonymous divergence)', size=14)
axs[0].set_ylabel('d$_N$/d$_S$', size=14, rotation='vertical')

# Set boundaries for the x and y axes of the main chart
axs[0].set_xlim([0.0000008, 0.1])
axs[0].set_ylim([0.05, 10])
axs[0].legend(fontsize=10)
axs[1].legend(fontsize=10,loc='upper right')

# Ensure a tight layout for the plot and show it
plt.tight_layout()
plt.show()
