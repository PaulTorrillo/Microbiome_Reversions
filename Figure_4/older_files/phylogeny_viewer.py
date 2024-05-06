import random
import Bio.Phylo
import numpy
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import SeqRecord, Seq
from Bio.Align import MultipleSeqAlignment
import re
import sys
import matplotlib.pyplot as plt
from Bio import Phylo
totalgens=5*int(20)#sys.argv[1])
nucleotides=['A','T','G','C']
stop_codons = ['TAA', 'TAG', 'TGA']
table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
codons = ['TTT', 'TTC', 'TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG',
          'ATT', 'ATC', 'ATA', 'ATG', 'GTT', 'GTC', 'GTA', 'GTG',
          'TCT', 'TCC', 'TCA', 'TCG', 'CCT', 'CCC', 'CCA', 'CCG',
          'ACT', 'ACC', 'ACA', 'ACG', 'GCT', 'GCC', 'GCA', 'GCG',
          'TAT', 'TAC','CAT', 'CAC', 'CAA', 'CAG',
          'AAT', 'AAC', 'AAA', 'AAG', 'GAT', 'GAC', 'GAA', 'GAG',
          'TGT', 'TGC', 'TGG', 'CGT', 'CGC', 'CGA', 'CGG',
          'AGT', 'AGC', 'AGA', 'AGG', 'GGT', 'GGC', 'GGA', 'GGG']


def classify_sites(sequence):
    num_sites = len(sequence)
    indices = list(range(num_sites))
    random.shuffle(indices)

    num_neutral = num_sites // 10
    num_beneficial = num_sites // 10

    neutral_sites = set(indices[:num_neutral])
    beneficial_sites = set(indices[num_neutral:num_neutral + num_beneficial])
    deleterious_sites = set(indices[num_neutral + num_beneficial:])

    return neutral_sites, beneficial_sites, deleterious_sites



def is_open_reading_frame(dna_seq):
    # Define start and stop codons
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']

    # Check if the sequence begins with a start codon
    if dna_seq.startswith(start_codon):
        # Find all possible stop codons
        possible_stops = [dna_seq[i:i + 3] for i in range(3, len(dna_seq)-3, 3)]

        # Check if there are no stop codons in between
        if all(codon not in stop_codons for codon in possible_stops):
            # Check if the sequence ends with a stop codon
            if dna_seq[-3:] in stop_codons:
                return True

    # If any of the checks fail, return False
    return False
def mutate(sequence,shadow_sequence, generations,numbs,neutral_sites, beneficial_sites, deleterious_sites):
    descendents=[]
    shadowdescendents = []
    if numbs>9:
        p=0
    else:
        p=min((10/numbs)**(1/(totalgens-generations))-1,1)
    sizeofsplit=random.choices([1,2], weights=[1-p,p], k=1)[0]
    for i in range(sizeofsplit):

        prob_skip_mutation=random.random()

        new_sequence = list(sequence)
        new_shadow_sequence = list(shadow_sequence)

        random_index = random.randint(0, len(new_sequence) - 1)
        codon = new_sequence[random_index]
        shadow_codon = new_shadow_sequence[random_index]

        # Choose a random nucleotide position in the codon to mutate
        nucleotide_position = random.randint(0, 2)
        original_nucleotide = codon[nucleotide_position]
        new_nucleotides = [n for n in ['A', 'T', 'G', 'C'] if n != original_nucleotide]
        new_nucleotide = random.choice(new_nucleotides)

        # Create new codons
        new_codon = codon[:nucleotide_position] + new_nucleotide + codon[nucleotide_position + 1:]
        new_shadow_codon = shadow_codon[:nucleotide_position] + new_nucleotide + shadow_codon[nucleotide_position + 1:]

        # Check if the new codon is a stop codon in the actual sequence
        if new_codon in stop_codons:
            descendents.append(new_sequence)
            shadowdescendents.append(new_shadow_sequence)
            continue # Skip the mutation if it results in a stop codon

        if new_shadow_codon in stop_codons:
            new_shadow_codon = shadow_codon  # Revert if it results in a stop codon

        # Check if mutation is synonymous
        original_aa = table[codon]
        new_aa = table[new_codon]

        if original_aa == new_aa and prob_skip_mutation<0.5:
            # Synonymous mutation
            new_sequence[random_index] = new_codon
            new_shadow_sequence[random_index] = new_shadow_codon
        else:
            # Non-synonymous mutation
            if random_index in neutral_sites and prob_skip_mutation<0.5:
                new_sequence[random_index] = new_codon
                new_shadow_sequence[random_index] = new_shadow_codon
            elif random_index in beneficial_sites:
                new_sequence[random_index] = new_codon
                if generations >= totalgens - 5:
                    new_shadow_sequence[random_index] = new_shadow_codon
                else:
                    # Find a codon in the shadow sequence that encodes the same amino acid
                    original_aa = table[shadow_codon]  # Get the amino acid of the original shadow codon
                    possible_nucleotides = ['A', 'T', 'G', 'C']
                    random.shuffle(possible_nucleotides)  # Shuffle the nucleotides

                    for nucleotide in possible_nucleotides:
                        potential_shadow_codon = (shadow_codon[:nucleotide_position] +
                                                  nucleotide +
                                                  shadow_codon[nucleotide_position + 1:])

                        # Check if the new codon encodes the same amino acid
                        if table[potential_shadow_codon] == original_aa:
                            new_shadow_sequence[random_index] = potential_shadow_codon
                            break
            elif random_index in deleterious_sites:
                # Attempt to find a non-synonymous mutation at a beneficial site
                mutation_found = False
                while not mutation_found and beneficial_sites:
                    beneficial_index = random.choice(list(beneficial_sites))
                    beneficial_codon = new_sequence[beneficial_index]
                    beneficial_shadow_codon = new_shadow_sequence[beneficial_index]
                    for attempt in range(10):
                        pos = random.randint(0, 2)
                        original_nuc = beneficial_codon[pos]
                        possible_nucs = [n for n in ['A', 'T', 'G', 'C'] if n != original_nuc]
                        new_nuc = random.choice(possible_nucs)
                        mutated_codon = beneficial_codon[:pos] + new_nuc + beneficial_codon[pos + 1:]
                        mutated_shadow_codon = beneficial_shadow_codon[:pos] + new_nuc + beneficial_shadow_codon[
                                                                                         pos + 1:]
                        if table[mutated_codon] != table[beneficial_codon] and mutated_codon not in stop_codons:
                            new_sequence[beneficial_index] = mutated_codon
                            if generations >= totalgens - 5:
                                if mutated_shadow_codon not in stop_codons:
                                    new_shadow_sequence[beneficial_index] = mutated_shadow_codon
                                else:
                                    new_shadow_sequence[beneficial_index] = beneficial_shadow_codon
                            else:
                                # Find a codon in the shadow sequence that encodes the same amino acid
                                original_aa = table[shadow_codon]  # Get the amino acid of the original shadow codon
                                possible_nucleotides = ['A', 'T', 'G', 'C']
                                random.shuffle(possible_nucleotides)  # Shuffle the nucleotides

                                for nucleotide in possible_nucleotides:
                                    potential_shadow_codon = (shadow_codon[:nucleotide_position] +
                                                              nucleotide +
                                                              shadow_codon[nucleotide_position + 1:])

                                    # Check if the new codon encodes the same amino acid
                                    if table[potential_shadow_codon] == original_aa:
                                        new_shadow_sequence[random_index] = potential_shadow_codon
                                        break
                            mutation_found = True
                            break
        descendents.append(new_sequence)
        shadowdescendents.append(new_shadow_sequence)
    return (descendents, shadowdescendents)

def generateRandomseqs(seqLength, generations):
    firstSequence=[]
    for i in range(seqLength):
        firstSequence=firstSequence+[random.choice(codons)]
    neutral_sites, beneficial_sites, deleterious_sites = classify_sites(firstSequence)
    currentleaves=[firstSequence]
    currentshadowleaves=[firstSequence]
    for j in range(generations):
        newleaves=[]
        shadowleaves=[]
        for i in range(len(currentleaves)):
            mutated_seq=mutate(currentleaves[i],currentshadowleaves[i], j, len(currentleaves),neutral_sites, beneficial_sites, deleterious_sites)
            list.extend(newleaves,mutated_seq[0])
            list.extend(shadowleaves,mutated_seq[1])
        currentleaves=newleaves
        currentshadowleaves=shadowleaves
    return (newleaves, shadowleaves)

zed=generateRandomseqs(500,totalgens)
seqs=[]
sseqs=[]
for i in range(len(zed[0])):
    seqs.append("ATG"+"".join(zed[0][i]))
    sseqs.append("ATG" + "".join(zed[1][i]))
# Align the sequences

records = [SeqRecord.SeqRecord(Seq.Seq(seq), id=f"seq_{i}") for i, seq in enumerate(seqs)]
srecords = [SeqRecord.SeqRecord(Seq.Seq(seq), id=f"seq_{i}") for i, seq in enumerate(sseqs)]
# Create a MultipleSeqAlignment object
alignment = MultipleSeqAlignment(records)
print(alignment)
salignment = MultipleSeqAlignment(srecords)
print(salignment)
# Print the alignment
with open('/Users/paultorrillo/Downloads/paml4.8/seqs.nuc', 'w') as f:
    print(alignment.__format__('phylip-sequential'), file=f)
with open('/Users/paultorrillo/Downloads/paml4.8/sseqs.nuc', 'w') as f:
    print(salignment.__format__('phylip-sequential'), file=f)
# Calculate the distance matrix
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(alignment)
sdm = calculator.get_distance(salignment)
# Construct the tree using the neighbor-joining algorithm
constructor = DistanceTreeConstructor(calculator, 'nj')
tree = constructor.build_tree(alignment)
tree.root_at_midpoint()

stree = constructor.build_tree(salignment)
stree.root_at_midpoint()

built_tree=tree.__format__('newick')
built_tree=built_tree.replace(":",' #')
built_tree=re.sub(r'Inner\d+', '', built_tree)
built_tree=built_tree.replace(" #0.00000",'')

shadow_built_tree=stree.__format__('newick')
shadow_built_tree=shadow_built_tree.replace(":",' #')
shadow_built_tree=re.sub(r'Inner\d+', '', shadow_built_tree)
shadow_built_tree=shadow_built_tree.replace(" #0.00000",'')

# Visualize and keep the "Normal Tree" open
plt.figure(1, figsize=(10, 5))  # Figure 1 for the normal tree
plt.title("Normal Tree")
Phylo.draw(tree)

# Visualize and keep the "Shadow Tree" open
plt.figure(2, figsize=(10, 5))  # Figure 2 for the shadow tree
plt.title("Shadow Tree")
Phylo.draw(stree)

# Use plt.show() to display both plots; this will keep the windows open
plt.show()