import numpy
import matplotlib
import matplotlib.pyplot as plt
import os
import scipy
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import ScalarFormatter
import pandas as pd
matplotlib.rcParams['font.family'] = 'Helvetica'
fig= plt.figure(figsize=(8, 8))

axs=plt.gca()

def to_year(x):
    years = x / (365)
    return years

def to_gen(x):
    gen=x*(365)
    return gen
# Define a function to format the ticks on the second axis
def years_formatter1(x, pos):
    years = int(x / 365)
    return f'{years}'

mutation_rate_per_bp_per_generation=10**-9
mutation_rate_per_codon_per_generation=3*mutation_rate_per_bp_per_generation
s=3.5*10**-5
core_genome_codons=500000
fraction_synonymous=0.25
fraction_nonsynonymous=0.75
fraction_N_neutral=1/10
neutral_N_rate_per_codon=mutation_rate_per_codon_per_generation*fraction_nonsynonymous*fraction_N_neutral
nonneutral_N_rate_per_codon=mutation_rate_per_codon_per_generation*fraction_nonsynonymous*(1-fraction_N_neutral)
S_rate_per_codon=mutation_rate_per_codon_per_generation*fraction_synonymous
normalization_constant=fraction_synonymous/fraction_nonsynonymous
numpoints=20000
max_generations=2000000
T_adapt=840
n_loci=55
#For top axis conversion

def to_year_2(x):
    years = x / (2 * S_rate_per_codon * 365)
    return years

def to_dS(x):
    dS=x*(2 * S_rate_per_codon * 365)
    return dS

def calculate_log_result(x, a, b):
    t = x / (2 * S_rate_per_codon)
    neutnon = a * mutation_rate_per_codon_per_generation * fraction_nonsynonymous * core_genome_codons * 2 * t
    syn = 2 * S_rate_per_codon * core_genome_codons * t
    t_mimic=numpy.full_like(t,10000)
    combine_t=[]
    combine_t.append(t)
    combine_t.append(t_mimic)
    newt=numpy.amin(combine_t,axis=0)
    result = (1 / 3) * ((neutnon + 2*b*newt)/ (syn))
    return numpy.log(result)
def small_round(x):
    counter=0
    newx=x
    while newx<100:
        newx=newx*10
        counter=counter+1
    newx=numpy.round(newx)
    newx=newx/(10**counter)
    return newx
def years_formatter2(x, pos):
    years = small_round((10**x)/(2*S_rate_per_codon*365))
    return f'{years:.1e}'
T_generations=numpy.linspace(1,max_generations,numpoints)
synonymous_muts=numpy.zeros(numpoints)
neutral_nonsynonymous_muts=numpy.zeros(numpoints)
dN_dS_theory=numpy.zeros(numpoints)
theory_fitness=numpy.zeros(numpoints)
Ns_sim_act=numpy.loadtxt("small_bottleneck_reversions.txt")
dN_dS_sim_obs=numpy.zeros(numpoints)
Ns_sim_obs=numpy.loadtxt("small_bottleneck_mutations.txt")
dN_dS_sim_act=numpy.zeros(numpoints)
Ns_hitch=numpy.loadtxt("small_bottleneck_hitch.txt")
dN_dS_1_quad=numpy.zeros(numpoints)
sim_fitness=numpy.ones(numpoints)

for n in range(0,numpoints):

    synonymous_muts[n]=S_rate_per_codon*core_genome_codons*2*T_generations[n]
    neutral_nonsynonymous_muts[n]=neutral_N_rate_per_codon*core_genome_codons*2*T_generations[n]
    dN_dS_theory[n]=normalization_constant*((((neutral_nonsynonymous_muts[n]+(n_loci*(1-numpy.exp((-2*T_generations[n])/(n_loci*T_adapt)))))/(synonymous_muts[n]))))
    theory_fitness[n] = (1-s)**(
                nonneutral_N_rate_per_codon * core_genome_codons * ((1 - numpy.exp(-s * T_generations[n])) / s))
    dN_dS_sim_obs[n] = normalization_constant * ((((neutral_nonsynonymous_muts[n] + 2 * Ns_sim_obs[n]) / (synonymous_muts[n]))))
    dN_dS_sim_act[n] = normalization_constant * ((((neutral_nonsynonymous_muts[n] + 2 * Ns_sim_act[n]) / (synonymous_muts[n]))))
    sim_fitness[n] = (1 - 0.003) ** Ns_hitch[n]


folder_path="../dnds_flat_files"
genuslist=set()
for file in os.listdir(folder_path):
    genuslist.add(file.split("_")[0]+'_'+file.split("_")[1])
genuslist=[*genuslist]
genuslist.sort()

for genus in genuslist:
    genustocombine=[]
    for file in os.listdir(folder_path):
        if file.startswith(genus):
            genustocombine.append(numpy.loadtxt('../dnds_flat_files/'+file,skiprows=1,delimiter=',',usecols=(0,1)))
    combined_genus=genustocombine[0]
    for i in range(len(genustocombine)):
        if i>0:
            combined_genus=numpy.concatenate((combined_genus, genustocombine[i]), axis=0)
    if numpy.sort(combined_genus[:, 0])[9]<0.0005:
        x=numpy.log(combined_genus[:, 0])
        y=numpy.log(combined_genus[:, 1])
        myline = numpy.linspace(numpy.log(0.000001), numpy.log(1), 1000)
        popt, pcov = scipy.optimize.curve_fit(calculate_log_result, numpy.exp(x), y)

axs.plot(T_generations, sim_fitness, linewidth=2, label='simulation')
axs.set_xlabel('bacterial generations',size=16)
axs.set_ylabel('relative fitness of population',size=16,rotation='vertical')
axs.tick_params(length=12, width=1,which='major',direction='inout', labelsize=14)
axs.tick_params(length=8, width=1,which='minor',direction='inout', labelsize=14)
axs.set_xlim([1,2000000])
# Configure the secondary x axis of the main chart
ax22 = axs.secondary_xaxis('top', functions=(to_year,to_gen))
ax22.set_xlabel('years', size=16)
ax22.tick_params(length=10, width=1, which='major', direction='inout',labelsize=14)
ax22.tick_params(length=6, width=1, which='minor', direction='inout', labelsize=14)

fig.tight_layout()
plt.show()
