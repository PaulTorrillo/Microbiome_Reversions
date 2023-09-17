# Import required libraries
import numpy
import matplotlib
import matplotlib.pyplot as plt
import os
import scipy
from matplotlib.ticker import FuncFormatter

# Set font for matplotlib
matplotlib.rcParams['font.family'] = 'Helvetica'

# Define mutation rates and proportions of synonymous and nonsynonymous codons
base_mutation_rate_per_generation=10**-9
codon_mutation_rate_per_generation=3*base_mutation_rate_per_generation
core_genome_codons=500000
synonymous_fraction=0.25
nonsynonymous_fraction=0.75
synonymous_rate_per_codon=codon_mutation_rate_per_generation*synonymous_fraction

# Define the logarithmic function for calculating the results
def log_result_func(x, s, a, b):
    intermediate_time = x / (2 * synonymous_rate_per_codon)
    neutnon = a * codon_mutation_rate_per_generation * nonsynonymous_fraction * core_genome_codons * 2 * intermediate_time
    syn = 2 * synonymous_rate_per_codon * core_genome_codons * intermediate_time
    result = (1 / 3) * ((neutnon + 2*b*((1-numpy.exp(-(2*intermediate_time)/s))/2))/ (syn))
    return numpy.log(result)

# Function for rounding to significant figures
def round_to_significant_figures(x):
    scaling_factor = 0
    scaled_x = x
    if scaled_x < 100:
        while scaled_x < 10:
            scaled_x *= 10
            scaling_factor += 1
        scaled_x = numpy.round(scaled_x)
        scaled_x /= 10**scaling_factor
    else:
        while scaled_x > 100:
            scaled_x /= 10
            scaling_factor += 1
        scaled_x = numpy.round(scaled_x)
        scaled_x *= 10**scaling_factor
    return scaled_x

# Function for formatting years in scientific notation
def format_years_scientifically(x, pos):
    years = round_to_significant_figures((10**x)/(2*synonymous_rate_per_codon*365))
    return f'{years:.1e}'

# Functions for performing optimization and returning specific parameters
def optimize_and_return_third_param(x,y):
    optimal_params, param_covariances = scipy.optimize.curve_fit(log_result_func, numpy.exp(x), y, bounds=([0, 0, 0], [numpy.inf, numpy.inf, numpy.inf]),p0=[1E5,0.1,20])
    return optimal_params[2]

def optimize_and_return_ratio_first_third_params(x,y):
    optimal_params, param_covariances = scipy.optimize.curve_fit(log_result_func, numpy.exp(x), y, bounds=([0, 0, 0], [numpy.inf, numpy.inf, numpy.inf]),p0=[1E5,0.1,20])
    return optimal_params[0]/optimal_params[2]

# Function to italicize text
def italicize(text):
    return r'$\it{%s}$' % text

# Directory containing the data files
data_folder_path="dnds_flat_files"

# Create a list of genus names from the file names
genus_list=set()
for filename in os.listdir(data_folder_path):
    genus_list.add(filename.split("_")[0] + '_' + filename.split("_")[1])
genus_list.remove('.DS_Store')  # Remove system file from list
genus_list = list(genus_list)
genus_list.sort()

# Create plot objects
fig, sub_plots = plt.subplots(1, 1, figsize=(8.5,4.2))
first_inset = sub_plots.inset_axes([0.67, 0.7, 0.32, 0.25])
second_inset = sub_plots.inset_axes([0.67, 0.4, 0.32, 0.25])

# Counter for tracking iteration
iteration_counter = 1

# Loop over the list of genera
for genus in genus_list:
    # Initialize a list to combine genus files
    genus_to_combine = []

    # Loop over the files in the folder_path
    for file in os.listdir(data_folder_path):
        # If the file starts with the current genus, append its data to the genus_to_combine list
        if file.startswith(genus):
            genus_to_combine.append(numpy.loadtxt(f'{data_folder_path}/{file}', skiprows=1, delimiter=',', usecols=(0,1)))

    # Combine all data for the current genus
    combined_genus_data = genus_to_combine[0]
    for i in range(1, len(genus_to_combine)):
        combined_genus_data = numpy.concatenate((combined_genus_data, genus_to_combine[i]), axis=0)

    # If the 10th smallest x-value is below a threshold, perform curve fitting and statistical tests
    if numpy.sort(combined_genus_data[:,0])[9] < 0.0005:
        x_values = numpy.log(combined_genus_data[:,0])
        y_values = numpy.log(combined_genus_data[:,1])

        line_space = numpy.linspace(numpy.log(0.000001), numpy.log(1), 1000)

        # Optimize curve fitting for the logarithmic function
        optimal_params, param_covariances = scipy.optimize.curve_fit(log_result_func, numpy.exp(x_values), y_values, p0=[1E5,0.1,20], bounds=([0, 0, 0], [numpy.inf, numpy.inf, numpy.inf]))

        # Calculate bootstrap estimates and plot confidence intervals
        bootstrap_result = scipy.stats.bootstrap((x_values, y_values), optimize_and_return_third_param, paired=True, n_resamples=999, method='basic')
        lower_confidence, upper_confidence = bootstrap_result.confidence_interval

        lower_error = min(optimal_params[2], optimal_params[2] - lower_confidence)
        upper_error = upper_confidence - optimal_params[2]

        second_inset.errorbar(iteration_counter, optimal_params[2], yerr=[[lower_error], [upper_error]], fmt='o', capsize=3)

        bootstrap_result2 = scipy.stats.bootstrap((x_values, y_values), optimize_and_return_ratio_first_third_params, paired=True, n_resamples=999, method='basic')
        lower_confidence2, upper_confidence2 = bootstrap_result2.confidence_interval

        lower_error2 = min((optimal_params[0]/optimal_params[2]) , (optimal_params[0]/optimal_params[2]) - lower_confidence2)
        upper_error2 = upper_confidence2 - (optimal_params[0]/optimal_params[2])

        first_inset.errorbar(iteration_counter, (optimal_params[0]/optimal_params[2]) , yerr=[[lower_error2], [upper_error2]], fmt='o', capsize=3)

        # Increment counter, scatter plot data and draw the fitted curve
        iteration_counter += 1
        sub_plots.scatter(numpy.exp(x_values), numpy.exp(y_values), s=10, alpha=0.25)
        sub_plots.plot(numpy.exp(line_space), numpy.exp(log_result_func(numpy.exp(line_space), *optimal_params)),
           label=italicize(genus).replace('Phocaeicola_', 'Ph. ').replace('Bacteroides_', 'B. ').replace('Parabacteroides_', 'Pa. '), linewidth='3')


# Set logarithmic scale for the main plot
sub_plots.set_xscale('log')
sub_plots.set_yscale('log')

# Set y-labels and limits for the inset plots
second_inset.set_ylabel('$n_{loci}$', rotation='vertical', size=10)
second_inset.set_ylim([-10, 110])

first_inset.set_ylim([-200, 2200])
first_inset.set_ylabel('$T_{adapt}$', rotation='vertical', size=10)

# Set tick parameters for the main plot
sub_plots.tick_params(length=10, width=1, which='major', direction='inout', labelsize=10)
sub_plots.tick_params(length=6, width=1, which='minor', direction='inout')

# Remove x-ticks from the inset plots and set their other tick parameters
second_inset.set_xticks([])
second_inset.tick_params(length=8, width=1, which='major', direction='inout', labelsize=10)
second_inset.tick_params(axis='x', which='both', bottom=False, top=False)

first_inset.set_xticks([])
first_inset.tick_params(length=8, width=1, which='major', direction='inout', labelsize=10)
first_inset.tick_params(axis='x', which='both', bottom=False, top=False)

# Set x- and y-labels for the main plot
sub_plots.set_xlabel('d$_S$ (core genome synonymous divergence)', size=12)
sub_plots.set_ylabel('d$_N$/d$_S$', size=12, rotation='vertical')

# Set x- and y-limits for the main plot
sub_plots.set_xlim([0.0000008, 0.1])
sub_plots.set_ylim([0.05, 10])

# Extract legend labels from the plot, and then set the legend location and other parameters
lines, labels = sub_plots.get_legend_handles_labels()
sub_plots.legend(lines, labels, loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=5, labelspacing=0.5, fontsize=10, frameon=False)

# Add a second x-axis to the top of the main plot, set its label, limit, formatter, and tick parameters
second_x_axis = sub_plots.twiny()
second_x_axis.set_xlabel('MRCA (years)', size=12, labelpad=10)
second_x_axis.set_xlim([-6.09, -1])
second_x_axis.xaxis.set_major_formatter(FuncFormatter(format_years_scientifically))
second_x_axis.tick_params(length=10, width=1, which='major', direction='inout', pad=12, labelsize=10)
second_x_axis.tick_params(length=6, width=1, which='minor', direction='inout')

# Adjust layout to make everything fit properly and then display the plot
plt.tight_layout()
plt.show()
