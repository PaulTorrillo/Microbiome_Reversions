import numpy

import matplotlib
import matplotlib.pyplot as plt

import os



from matplotlib.ticker import FuncFormatter
from matplotlib.gridspec import GridSpec

matplotlib.rcParams['font.family'] = 'Helvetica'

fig = plt.figure(figsize=(8.5, 8))

gs = GridSpec(2, 1, figure=fig, height_ratios=[1, 2])

# Add the first subplot on the top using the GridSpec gs[0]
ax2 = fig.add_subplot(gs[0])

# Add the second subplot at the bottom using the GridSpec gs[1:]
ax3 = fig.add_subplot(gs[1:])

num_points = 10000  # 5 years worth of data
x_days = numpy.arange(num_points) * 10  # days


x_years = x_days / 365

# Create a second axis for years



# Define a function to format the ticks on the second axis
def years_formatter1(x, pos):
    years = int(x / 365)
    return f'{years}'

mutation_rate_per_bp_per_generation=10**-9

mutation_rate_per_codon_per_generation=3*mutation_rate_per_bp_per_generation

s=3.5*10**-5

codons_per_genome=10**6

core_genome_codons=500000#0.4*codons_per_genome

fraction_synonymous=0.25

fraction_nonsynonymous=0.75

fraction_N_neutral=1/10

neutral_N_rate_per_codon=mutation_rate_per_codon_per_generation*fraction_nonsynonymous*fraction_N_neutral

nonneutral_N_rate_per_codon=mutation_rate_per_codon_per_generation*fraction_nonsynonymous*(1-fraction_N_neutral)

S_rate_per_codon=mutation_rate_per_codon_per_generation*fraction_synonymous

print(neutral_N_rate_per_codon*500000)

normalization_constant=fraction_synonymous/fraction_nonsynonymous

numpoints=50000

max_generations=5000000


def small_round(x):
    counter=0
    newx=x
    while newx<100:
        newx=newx*10
        counter=counter+1
    newx=numpy.round(newx)
    newx=newx/(10**counter)
    return newx

def years_formatter(x, pos):
    years = small_round((10**x)/(2*S_rate_per_codon*365))
    return f'{years:.1e}'

T_generations=numpy.linspace(1,max_generations,numpoints)

synonymous_muts=numpy.zeros(numpoints)

neutral_nonsynonymous_muts=numpy.zeros(numpoints)

dN_dS_theory=numpy.zeros(numpoints)
theory_fitness=numpy.zeros(numpoints)


Ns_sim= numpy.loadtxt("ultramutation2.txt")
dN_dS_sim=numpy.zeros(numpoints)
sim_fitness=numpy.ones(numpoints)


for n in range(0,numpoints):

    synonymous_muts[n]=S_rate_per_codon*core_genome_codons*2*T_generations[n]


    neutral_nonsynonymous_muts[n]=neutral_N_rate_per_codon*core_genome_codons*2*T_generations[n]


    dN_dS_theory[n]=normalization_constant*((((neutral_nonsynonymous_muts[n]+2*nonneutral_N_rate_per_codon*core_genome_codons*((1-numpy.exp(-s*T_generations[n]))/s))/(synonymous_muts[n]))))
    theory_fitness[n] = (1-s)**(
                nonneutral_N_rate_per_codon * core_genome_codons * ((1 - numpy.exp(-s * T_generations[n])) / s))

    dN_dS_sim[n] = normalization_constant * ((((2 * Ns_sim[n]) / (synonymous_muts[n]))))

    sim_fitness[n] = (1 - s) ** Ns_sim[n]


ax3.plot(S_rate_per_codon*2*T_generations,dN_dS_theory,linewidth=3,label='theory',color='k')


ax3.plot(S_rate_per_codon * 2 * T_generations, dN_dS_sim, linewidth=3, label='simulation')


ax3.plot(S_rate_per_codon*2*T_generations,numpy.ones_like(dN_dS_theory),linewidth=3,label='d$_N$/d$_S$ = 1',linestyle='--',color='k')

ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.tick_params(length=12, width=1,which='major',direction='inout')
ax3.tick_params(length=8, width=1,which='minor',direction='inout')
ax3.set_xlabel('d$_S$ (core genome synonymous divergence)',size=14)
ax3.set_ylabel('d$_N$/d$_S$',size=14,rotation='vertical')
ax3.set_xlim([0.00001,0.01])
ax3.set_ylim([0.05,5])
ax3.legend(fontsize='10',ncol=2,frameon=False)



folder_path="dnds_flat_files"
genuslist=set()
for file in os.listdir(folder_path):
    genuslist.add(file.split("_")[0]+'_'+file.split("_")[1])
print(genuslist)
genuslist.remove('.DS_Store')
genuslist=[*genuslist]
genuslist.sort()


for genus in genuslist:
    genustocombine=[]
    for file in os.listdir(folder_path):
        if file.startswith(genus):
            genustocombine.append(numpy.loadtxt('dnds_flat_files/'+file,skiprows=1,delimiter=',',usecols=(0,1)))
    combined_genus=genustocombine[0]
    for i in range(len(genustocombine)):
        if i>0:
            combined_genus=numpy.concatenate((combined_genus, genustocombine[i]), axis=0)

    a= numpy.argsort(combined_genus[:, 1])[::-1][0:10]

    b = [i for i, v in enumerate(combined_genus[:, 0]) if v < 0.0006]
    if numpy.sort(combined_genus[:, 0])[9]<0.0005:
        x=numpy.log(combined_genus[:, 0])
        y=numpy.log(combined_genus[:, 1])
        ax3.scatter(numpy.exp(x), numpy.exp(y),s=5,alpha=0.2,color='gray')

ax32=ax3.twiny()
ax32.set_xlabel('MRCA (years)',size=14,labelpad=10)


# Define a function to format the ticks on the second axis


    # Set the tick formatter for the second axis
ax32.set_xlim([-5,-2])
ax32.xaxis.set_major_formatter(FuncFormatter(years_formatter))
ax32.tick_params(length=12, width=1,which='major',direction='inout')
ax32.tick_params(length=8, width=1,which='minor',direction='inout')


ax2.plot(T_generations,theory_fitness,linewidth=3,label='theory',color='k')


ax2.plot(T_generations, sim_fitness, linewidth=3, label='simulation')





ax2.set_xlabel('bacterial generations',size=14)
ax2.set_ylabel('relative fitness',size=14,rotation='vertical')
ax2.tick_params(length=12, width=1,which='major',direction='inout')
ax2.tick_params(length=8, width=1,which='minor',direction='inout')
ax2.legend(fontsize=10,frameon=False)
ax2.set_xlim([1,500000])
ax2.set_ylim([0.99,1])
ax22 = ax2.twiny()

# Set the limits of the second axis in years
ax22.set_xlim(ax2.get_xlim())
ax22.set_xlabel('years', size=14,labelpad=10)




# Set the tick formatter for the second axis
ax22.xaxis.set_major_formatter(FuncFormatter(years_formatter1))

# Show the plot
# matplotlib.pyplot.yscale('log')
ax22.tick_params(length=12, width=1, which='major', direction='inout')
ax22.tick_params(length=8, width=1, which='minor', direction='inout')

fig.tight_layout()

plt.show()