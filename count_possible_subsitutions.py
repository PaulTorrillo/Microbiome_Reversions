from Bio.Data import CodonTable

def calculate_mutations(codon_table):
    bases = ['A', 'C', 'G', 'T']
    synonymous = 0
    nonsynonymous = 0
    stop_codon_mutations = 0  # Counter for mutations leading to stop codons

    for codon in codon_table.forward_table.keys():  # Iterate through all codons that encode an amino acid
        amino_acid = codon_table.forward_table[codon]
        for i in range(3):  # For each position in the codon
            for base in bases:  # For each possible mutation
                if codon[i] != base:  # If it's a mutation
                    mutated_codon = codon[:i] + base + codon[i + 1:]  # Create the mutated codon
                    # Determine if the mutation is synonymous, nonsynonymous, or leads to a stop codon
                    if mutated_codon in codon_table.forward_table:
                        mutated_amino_acid = codon_table.forward_table[mutated_codon]
                        if amino_acid == mutated_amino_acid:
                            synonymous += 1
                        else:
                            nonsynonymous += 1
                    elif mutated_codon in codon_table.stop_codons:  # Check if it's a stop codon
                        stop_codon_mutations += 1
                    else:
                        nonsynonymous += 1  # In case it's an unknown mutation

    return synonymous, nonsynonymous, stop_codon_mutations

# Using the bacterial codon table (Table 11)
bacterial_codon_table = CodonTable.unambiguous_dna_by_id[11]

synonymous_mutations, nonsynonymous_mutations, stop_codon_mutations = calculate_mutations(bacterial_codon_table)

print('Synonymous: '+str(synonymous_mutations))
print('Nonsynonymous excluding stop: '+str(nonsynonymous_mutations))
print('Stop codon: '+str(stop_codon_mutations))
