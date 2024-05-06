# Import required libraries
import numpy
import matplotlib
import matplotlib.pyplot as plt
import os
import scipy
from matplotlib.ticker import FuncFormatter
from matplotlib.gridspec import GridSpec
# Set font for matplotlib
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
def calculate_result(x, s, a, b):
    intermediate_time = x / (2 * synonymous_rate_per_codon)
    neutnon = a * codon_mutation_rate_per_generation * nonsynonymous_fraction * core_genome_codons * 2 * intermediate_time
    syn = 2 * synonymous_rate_per_codon * core_genome_codons * intermediate_time
    result = (1 / 3) * ((neutnon + 2*b*((1-numpy.exp(-(2*intermediate_time)/s))/2))/ (syn))
    return result

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
    optimal_params, param_covariances = scipy.optimize.curve_fit(calculate_result, x, y, bounds=([0, 0, 0], [numpy.inf, numpy.inf, numpy.inf]), p0=[1E5, 0.1, 20])
    return optimal_params[2]

def optimize_and_return_ratio_first_third_params(x,y):
    optimal_params, param_covariances = scipy.optimize.curve_fit(calculate_result, x, y, bounds=([0, 0, 0], [numpy.inf, numpy.inf, numpy.inf]), p0=[1E5, 0.1, 20])
    return optimal_params[0]

# Function to italicize text
def italicize(text):
    return r'$\it{%s}$' % text

# Directory containing the data files
data_folder_path="../dnds_flat_files"
collected_rs=[]

# Create a list of genus names from the file names
genus_list=set()
for filename in os.listdir(data_folder_path):
    genus_list.add(filename.split("_")[0] + '_' + filename.split("_")[1])
genus_list = list(genus_list)
genus_list.sort()

# Create plot objects
# Plot settings
fig, axs = plt.subplots(3, 1, figsize=(8.5, 8.2), gridspec_kw={'height_ratios': [1.9,  0.55,0.55]})
ax_inset = axs[0].inset_axes([0.66, 0.54, 0.3, 0.4])

# Counter for tracking iteration
iteration_counter = 1
findingmedianT=[]
findingmediann=[]
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
        x_values = combined_genus_data[:,0]
        y_values = combined_genus_data[:,1]

        line_space = numpy.linspace(numpy.log(0.000001), numpy.log(1), 1000)

        # Optimize curve fitting for the logarithmic function
        optimal_params, param_covariances = scipy.optimize.curve_fit(calculate_result, x_values, y_values, p0=[1E5, 0.1, 20], bounds=([0, 0, 0], [numpy.inf, numpy.inf, numpy.inf]))

        # Calculate bootstrap estimates and plot confidence intervals
        bootstrap_result = scipy.stats.bootstrap((x_values, y_values), optimize_and_return_third_param, paired=True, n_resamples=999, method='basic')
        lower_confidence, upper_confidence = bootstrap_result.confidence_interval

        lower_error = min(optimal_params[2], optimal_params[2] - lower_confidence)
        upper_error = upper_confidence - optimal_params[2]

        axs[2].errorbar(iteration_counter, optimal_params[2], yerr=[[lower_error], [upper_error]], fmt='o', capsize=3)

        findingmediann.append(optimal_params[2])
        print(optimal_params[2])
        bootstrap_result2 = scipy.stats.bootstrap((x_values, y_values), optimize_and_return_ratio_first_third_params, paired=True, n_resamples=999, method='basic')
        lower_confidence2, upper_confidence2 = bootstrap_result2.confidence_interval

        lower_error2 = min((optimal_params[0]) , (optimal_params[0]) - lower_confidence2)
        upper_error2 = upper_confidence2 - (optimal_params[0])

        axs[1].errorbar(iteration_counter, (optimal_params[0]) , yerr=[[lower_error2], [upper_error2]], fmt='o', capsize=3)
        findingmedianT.append(optimal_params[0])
        # Increment counter, scatter plot data and draw the fitted curve
        iteration_counter += 1
        axs[0].scatter(x_values, y_values, s=10, alpha=0.25)
        ax_inset.scatter(x_values, y_values * x_values, s=10, alpha=0.25)
        residuals = y_values - calculate_result(x_values, *optimal_params)
        var = y_values - numpy.mean(y_values)
        rss = numpy.sum(residuals ** 2)
        tss = numpy.sum(var ** 2)
        z = round(1 - (rss / tss), 3)
        collected_rs.append(z)
        axs[0].plot(numpy.exp(line_space), calculate_result(numpy.exp(line_space), *optimal_params),
                    label=italicize(genus).replace('Phocaeicola_', 'Ph. ').replace('Bacteroides_', 'B. ').replace(
                        'Parabacteroides_', 'Pa. '),
                    linewidth='3')

        ax_inset.plot(numpy.exp(line_space), calculate_result(numpy.exp(line_space), *optimal_params) * numpy.exp(line_space))

# Set logarithmic scale for the main plot
axs[0].set_xscale('log')
axs[0].set_yscale('log')
# Set boundaries for the inset chart
ax_inset.set_xlim([0.00001, 0.0003])
# Set boundaries for the y axis of the inset chart
ax_inset.set_ylim([0.00001, 0.0001])

# Set labels for the inset chart axes
ax_inset.set_xlabel('d$_S$', size=12)
ax_inset.set_ylabel('d$_N$', size=12)

# Set tick parameters for the inset and main charts
ax_inset.tick_params(length=8, width=1, which='major', direction='inout', labelsize=10)
# Set y-labels and limits for the inset plots
axs[2].set_ylabel('n$_{loci}$', rotation='vertical', size=10)
#axs[2].set_ylim([-15, 165])


axs[1].set_ylabel(r'$\tau_{flip}$', rotation='vertical', size=10)

# Set tick parameters for the main plot
axs[0].tick_params(length=10, width=1, which='major', direction='inout', labelsize=10)
axs[0].tick_params(length=6, width=1, which='minor', direction='inout')

# Remove x-ticks from the inset plots and set their other tick parameters
axs[2].set_xticks([])
axs[2].tick_params(length=8, width=1, which='major', direction='inout', labelsize=10)
axs[2].tick_params(axis='x', which='both', bottom=False, top=False)

axs[1].set_xticks([])
axs[1].tick_params(length=8, width=1, which='major', direction='inout', labelsize=10)
axs[1].tick_params(axis='x', which='both', bottom=False, top=False)

# Set x- and y-labels for the main plot
axs[0].set_xlabel('d$_S$ (core genome synonymous divergence)', size=12)
axs[0].set_ylabel('d$_N$/d$_S$', size=12, rotation='vertical')

# Set x- and y-limits for the main plot
axs[0].set_xlim([0.0000008, 0.1])
axs[0].set_ylim([0.05, 10])

# Extract legend labels from the plot, and then set the legend location and other parameters
lines, labels = axs[0].get_legend_handles_labels()
axs[0].legend(lines, labels, loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=5, labelspacing=0.5, fontsize=10, frameon=False)

# Configure the secondary x axis of the main chart
ax2 = axs[0].secondary_xaxis('top', functions=(to_year,to_dS))
ax2.set_xlabel('MRCA (years)', size=14)
ax2.tick_params(length=10, width=1, which='major', direction='inout',labelsize=12)
ax2.tick_params(length=6, width=1, which='minor', direction='inout')

print(numpy.max(collected_rs))
print(numpy.min(collected_rs))
# Adjust layout to make everything fit properly and then display the plot
plt.tight_layout()
plt.show()
