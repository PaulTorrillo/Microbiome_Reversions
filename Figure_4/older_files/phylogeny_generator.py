import random
import Bio.Phylo
import numpy
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import SeqRecord, Seq
from Bio.Align import MultipleSeqAlignment
import re
import sys
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
def mutate(sequence,shadowsequence, generations,numbs):
    descendents=[]
    shadowdescendents = []
    if numbs>9:
        p=0
    else:
        p=min((10/numbs)**(1/(totalgens-generations))-1,1)
    sizeofsplit=random.choices([1,2], weights=[1-p,p], k=1)[0]
    for i in range(sizeofsplit):
        newsequence=list(sequence)
        newshadowsequence = list(shadowsequence)
        random_index1 = random.randint(0, len(newsequence)-1)
        random_index2 = random.randint(0, 2)
        nucs=['A','T','G','C']
        invalidmutation=True
        currentcodon=newsequence[random_index1]
        nucs.remove(list(currentcodon)[random_index2])
        changing=''
        while invalidmutation:
            potetialswitch = list(currentcodon)
            changing=random.choice(nucs)
            potetialswitch[random_index2]=changing
            potetialswitch="".join(potetialswitch)
            if potetialswitch not in stop_codons:
                invalidmutation=False
        aa1=table[currentcodon]
        aa2=table[potetialswitch]
        rando2=random.uniform(0,1)
        if (aa1==aa2):
            if rando2<(1/3):
                newsequence[random_index1]=potetialswitch
                newshadowsequence[random_index1] = potetialswitch
        elif ((random.uniform(0, 1)<0.1)):
            if rando2<1/3:
                newsequence[random_index1] = potetialswitch
                newshadowsequence[random_index1] = potetialswitch
        else:
            newshadowsequence[random_index1] = potetialswitch
            if generations>=totalgens-5:
                newsequence[random_index1] = potetialswitch
        descendents.append(newsequence)
        shadowdescendents.append(newshadowsequence)
    return (descendents,shadowdescendents)



def generateRandomseqs(seqLength, generations):
    firstSequence=[]
    for i in range(seqLength):
        firstSequence=firstSequence+[random.choice(codons)]
    currentleaves=[firstSequence]
    currentshadowleaves=[firstSequence]
    for j in range(generations):
        newleaves=[]
        shadowleaves=[]
        for i in range(len(currentleaves)):
            mutated_seq=mutate(currentleaves[i],currentshadowleaves[i], j, len(currentleaves))
            list.extend(newleaves,mutated_seq[0])
            list.extend(shadowleaves,mutated_seq[1])
        currentleaves=newleaves
        currentshadowleaves=shadowleaves
    return (newleaves, shadowleaves)

zed=generateRandomseqs(330,totalgens)
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
salignment = MultipleSeqAlignment(srecords)
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
# Draw the tree using a graphical interface
with open('/Users/paultorrillo/Downloads/paml4.8/normal.treefile', 'w') as f:
    print(built_tree, file=f)
shadow_built_tree=stree.__format__('newick')
shadow_built_tree=shadow_built_tree.replace(":",' #')
shadow_built_tree=re.sub(r'Inner\d+', '', shadow_built_tree)
shadow_built_tree=shadow_built_tree.replace(" #0.00000",'')
# Draw the tree using a graphical interface
with open('/Users/paultorrillo/Downloads/paml4.8/shadow.treefile', 'w') as f:
    print(shadow_built_tree, file=f)
