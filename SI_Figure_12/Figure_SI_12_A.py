import numpy as np
import matplotlib.pyplot as plt
import os
import scipy
from matplotlib.ticker import FuncFormatter
import matplotlib
from scipy.sparse import lil_matrix  # LIL format is good for incremental construction
from scipy.sparse.linalg import spsolve  # For solving sparse linear systems
from scipy import sparse
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
n_loci=55
T_adapt=840
time_step = 100/(n_loci*T_adapt)
number = 99999
n = int(2 * time_step * number) + 5

def to_year(x):
    years = x / (2 * synonymous_rate_per_codon * 365)
    return years

def to_dS(x):
    dS=x*(2 * synonymous_rate_per_codon * 365)
    return dS

# Function for calculating log re
def calculate_weighted_sum(p):
    # Create an n x n sparse matrix
    A = lil_matrix((n, n))

    # Construct the matrix based on the described rules
    A[1, 0] = time_step  # First row transitions to the second row with probability 1
    A[0, 0] = 1 - time_step
    for i in range(1, n - 1):  # Apply rules for intermediate rows
        A[i - 1, i] = (1 - p) * time_step  # Backward transition
        A[i + 1, i] = p * time_step  # Forward transition
        A[i, i] = 1 - time_step
    A[n - 2, n - 1] = time_step  # Last row transitions backwards to the second-to-last row with probability 1
    A[n - 2, n - 2] = 1 - time_step  # 0.5

    # Convert A to CSR format for faster arithmetic operations
    A = A.tocsr()

    # Initial vector with 1 in the first row and 0 in all other rows
    v = np.array([1 if i == 0 else 0 for i in range(n)])

    # Multiply the vector by the matrix 100 times
    resultant_list = []
    indices = np.arange(n)
    for i in range(number):
        v = A.dot(v)
        v = v / np.sum(v)
        resultant_list.append((n_loci*np.dot(indices, v) / (3*(i + 1)*100*500000*7.5*10**-10)) + 0.1)

    weighted_sum = np.dot(indices, v)

    return resultant_list

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
    years = small_round((10 ** x) / (synonymous_rate_per_codon * 365))
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
#genus_list.remove('.DS_Store')
genus_list = [*genus_list]
genus_list.sort()

# Plot settings
fig, axs = plt.subplots(1, 1, figsize=(8.5, 4.1))#, gridspec_kw={'height_ratios': [1, 1]})

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

        genus_counter = genus_counter + 1
        axs.scatter(np.exp(x), np.exp(y), s=5, alpha=0.20,color='gray')
        #axs[1].scatter(np.exp(x), np.exp(y) * np.exp(x), s=10, alpha=0.25,color='gray')


indices = np.arange(number)*200*7.5*10**-10
plot_x=calculate_weighted_sum(0)
axs.plot(indices, plot_x,
               label='0% compensatory, n$_{loci}=55$',
              linewidth='3')
#axs[1].plot(indices, plot_x*indices,
 #              label='100% true reversion, n$_{loci}=55$',
  #            linewidth='3')

indices = np.arange(number)*200*7.5*10**-10
plot_x=calculate_weighted_sum(0.25)
axs.plot(indices, plot_x,
               label='25% compensatory, n$_{loci}=55$',
              linewidth='3')
#axs[1].plot(indices, plot_x*indices,
 #              label='75% true reversion, n$_{loci}=55$',
  #            linewidth='3')
indices = np.arange(number)*200*7.5*10**-10
plot_x=calculate_weighted_sum(0.5)
axs.plot(indices, plot_x,
               label='50% compensatory, n$_{loci}=55$',
              linewidth='3')
#axs[1].plot(indices, plot_x*indices,
 #              label='50% true reversion, n$_{loci}=55$',
        #      linewidth='3')
indices = np.arange(number)*200*7.5*10**-10
plot_x=calculate_weighted_sum(0.75)
axs.plot(indices, plot_x,
               label='75% compensatory, n$_{loci}=55$',
              linewidth='3')
#axs[1].plot(indices, plot_x*indices,
  #             label='25% true reversion, n$_{loci}=55$',
  #            linewidth='3')
indices = np.arange(number)*200*7.5*10**-10
n_loci=28
time_step = 100/(n_loci*T_adapt)
plot_x=calculate_weighted_sum(0.25)
axs.plot(indices, plot_x,
               label='25% compensatory, n$_{loci}=28$',
              linewidth='3',linestyle='--',color='tab:orange')

 #              label='75% true reversion, n$_{loci}=28$',
  #            linewidth='3',linestyle='--',color='tab:orange')
# Set boundaries for the inset chart
#axs[1].set_xlim([0.0000, 0.0006])

# Set axes of main plot to logarithmic scale
axs.set_xscale('log')
axs.set_yscale('log')



# Set boundaries for the y axis of the inset chart
#axs[1].set_ylim([0.0000, 0.0002])

# Set labels for the inset chart axes
#axs[1].set_xlabel('d$_S$', size=12)
#axs[1].set_ylabel('d$_N$', size=12)

# Set tick parameters for the inset and main charts
#xs[1].tick_params(length=8, width=1, which='major', direction='inout', labelsize=10)
axs.tick_params(length=10, width=1, which='major', direction='inout', labelsize=12)
axs.tick_params(length=6, width=1, which='minor', direction='inout')


# Set labels for the x and y axes of the main chart
axs.set_xlabel('d$_S$ (core genome synonymous divergence)', size=14)
axs.set_ylabel('d$_N$/d$_S$', size=14, rotation='vertical')

# Set boundaries for the x and y axes of the main chart
axs.set_xlim([0.0000008, 0.01])
axs.set_ylim([0.05, 10])
axs.legend(fontsize=8,loc='upper right')
#axs[1].legend(fontsize=10,loc='lower right')
# Configure the secondary x axis of the main chart
ax2 = axs.secondary_xaxis('top', functions=(to_year,to_dS))
ax2.set_xlabel('MRCA (years)', size=14)
ax2.tick_params(length=10, width=1, which='major', direction='inout',labelsize=12)
ax2.tick_params(length=6, width=1, which='minor', direction='inout')
# Ensure a tight layout for the plot and show it
plt.tight_layout()
plt.show()