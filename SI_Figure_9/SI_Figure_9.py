import numpy as np
import matplotlib.pyplot as plt
import os
import scipy
from matplotlib.ticker import FuncFormatter
import matplotlib

# Configuring matplotlib settings
matplotlib.rcParams['font.family'] = 'Helvetica'

# Define mutation rates and proportions of synonymous and nonsynonymous codons
base_mutation_rate_per_generation=10**-9
codon_mutation_rate_per_generation=3*base_mutation_rate_per_generation
core_genome_codons=500000
synonymous_fraction=0.25
nonsynonymous_fraction=0.75
synonymous_rate_per_codon=codon_mutation_rate_per_generation*synonymous_fraction
def to_year(x):
    years = x / (2 * synonymous_rate_per_codon * 365)
    return years

def to_dS(x):
    dS=x*(2 * synonymous_rate_per_codon * 365)
    return dS

# Define the logarithmic function for calculating the results
def log_result_func(x, s, a, b):
    intermediate_time = x / (2 * synonymous_rate_per_codon)
    neutnon = a * codon_mutation_rate_per_generation * nonsynonymous_fraction * core_genome_codons * 2 * intermediate_time
    syn = 2 * synonymous_rate_per_codon * core_genome_codons * intermediate_time
    result = (1 / 3) * ((neutnon + 2*b*((1-np.exp(-(2*intermediate_time)/s))/2))/ (syn))
    return np.log(result)

# Function for rounding to significant figures
def round_to_significant_figures(x):
    scaling_factor = 0
    scaled_x = x
    if scaled_x < 100:
        while scaled_x < 10:
            scaled_x *= 10
            scaling_factor += 1
        scaled_x = np.round(scaled_x)
        scaled_x /= 10**scaling_factor
    else:
        while scaled_x > 100:
            scaled_x /= 10
            scaling_factor += 1
        scaled_x = np.round(scaled_x)
        scaled_x *= 10**scaling_factor
    return scaled_x

# Function for formatting years in scientific notation
def format_years_scientifically(x, pos):
    years = round_to_significant_figures((10**x)/(2*synonymous_rate_per_codon*365))
    return f'{years:.1e}'

# Functions for performing optimization and returning specific parameters
def optimize_and_return_third_param(x,y):
    optimal_params, param_covariances = scipy.optimize.curve_fit(log_result_func, np.exp(x), y, bounds=([0, 0, 0], [np.inf, np.inf, np.inf]),p0=[1E5,0.1,20])
    return optimal_params[2]

def optimize_and_return_ratio_first_third_params(x,y):
    optimal_params, param_covariances = scipy.optimize.curve_fit(log_result_func, np.exp(x), y, bounds=([0, 0, 0], [np.inf, np.inf, np.inf]),p0=[1E5,0.1,20])
    return optimal_params[0]

# Function to italicize text
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
    if np.sort(combined_genus[:, 0])[10] < 0.0001:
        x_values = np.log(combined_genus[:, 0])
        y_values = np.log(combined_genus[:, 1])

        line_space = np.linspace(np.log(0.000001), np.log(1), 1000)

        # Optimize curve fitting for the logarithmic function
        optimal_params, param_covariances = scipy.optimize.curve_fit(log_result_func, np.exp(x_values), y_values,
                                                                     p0=[1E5, 0.1, 20], bounds=(
            [0, 0, 0], [np.inf, np.inf, np.inf]))


        # Increment counter, scatter plot data and draw the fitted curve
        genus_counter += 1
        color = cmap(genus_counter)
        plt.scatter(np.exp(x_values), np.exp(y_values), s=10, alpha=0.25,color=color)
        plt.plot(np.exp(line_space), np.exp(log_result_func(np.exp(line_space), *optimal_params)),
                       label=italicize(genus).replace('Phocaeicola_', 'Ph. ').replace('Bacteroides_', 'B. ').replace(
                           'Parabacteroides_', 'Pa. '), linewidth='3',color=color)

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
