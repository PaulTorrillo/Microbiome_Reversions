import numpy as np
import scipy.stats
from Bio import SeqIO
import scipy.stats as stats
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import ks_2samp
from scipy.stats import ttest_ind


cog_categories = {
    'J': 'Translation, ribosomal structure and biogenesis',
    'A': 'RNA processing and modification',
    'K': 'Transcription',
    'L': 'Replication, recombination and repair',
    'B': 'Chromatin structure and dynamics',
    'D': 'Cell cycle control, cell division, chromosome partitioning',
    'Y': 'Nuclear structure',
    'V': 'Defense mechanisms',
    'T': 'Signal transduction mechanisms',
    'M': 'Cell wall/membrane/envelope biogenesis',
    'N': 'Cell motility',
    'Z': 'Cytoskeleton',
    'W': 'Extracellular structures',
    'U': 'Intracellular trafficking, secretion, and vesicular transport',
    'O': 'Posttranslational modification, protein turnover, chaperones',
    'C': 'Energy production and conversion',
    'G': 'Carbohydrate transport and metabolism',
    'E': 'Amino acid transport and metabolism',
    'F': 'Nucleotide transport and metabolism',
    'H': 'Coenzyme transport and metabolism',
    'I': 'Lipid transport and metabolism',
    'P': 'Inorganic ion transport and metabolism',
    'Q': 'Secondary metabolites biosynthesis, transport, and catabolism',
    'R': 'General function prediction only',
    'S': 'Function unknown'
}


def map_characters_to_entries(input_string, char_to_entry_dict):
    # Initialize an empty list to hold the mapped entries
    mapped_entries = []

    # Iterate over each character in the input string
    for char in input_string:
        # Use the dictionary to get the corresponding entry for the current character
        # If the character is not in the dictionary, you can choose to skip it or handle it differently
        entry = char_to_entry_dict.get(char)

        # If an entry is found, append it to the list of mapped entries
        if entry is not None:
            mapped_entries.append(entry)

    # Join the list of mapped entries with '/' and return the result
    return '/'.join(mapped_entries)
def find_string_in_excel(file_path, search_str):
    # Load the Excel file
    df = pd.read_excel(file_path)
    print(search_str)
    # Check if the search string is in the first column
    result_row = df[df[df.columns[0]] == search_str]

    if not result_row.empty:
        # Assuming the second column is what you want to return
        return result_row.iloc[0, 1]
    else:
        return "String not found in the first column."
def remove_first_word_preserve_whitespace(s):
    # Find the index of the first space
    first_space_index = s.find(' ')

    # If there's no space, return the whole string
    if first_space_index == -1:
        return ''

    # Return the substring from the first space onwards, preserving original whitespace
    return s[first_space_index + 1:]

def calculate_gc_content(sequence):
    """Calculate the GC content of a DNA sequence."""
    # Ensure the sequence is uppercase for consistency
    sequence = sequence.upper()
    # Count the G and C nucleotides in the sequence
    gc_count = sequence.count('G') + sequence.count('C')
    # Calculate the GC content percentage
    gc_content = (gc_count / len(sequence))
    return gc_content
def smooth_data_with_level(data, level=5):
    """
    Smooths the data based on a given level of smoothing.
    The level determines how many elements to the left and right are considered for averaging.

    :param data: List of numbers to be smoothed.
    :param level: Number of elements to the left and right to include in the averaging.
    :return: List of smoothed numbers.
    """
    smoothed_data = []
    for i in range(len(data)):
        # Calculate the indices of the neighboring elements
        left = max(0, i - level)
        right = min(len(data), i + level + 1)  # +1 because the range's stop is exclusive
        # Calculate the average
        average = sum(data[left:right]) / (right - left)
        smoothed_data.append(average)
    return smoothed_data



def run_main_code(file_path):
    leucine_codons = {'TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'}
    serine_codons = {'TCA', 'TCT', 'TCG', 'TCC', 'AGC', 'AGT'}
    special_codons = {'TTA', 'TCA'}  # TTA for Leucine, TCA for Serine

    gene_proportions = {}
    susC_gene_proportions = {}
    plotting_loci=[]
    plotting_gc=[]
    plotting_deviations=[]
    plotting_name=[]
    plotting_size=[]
    plotting_location=[]
    plotting_start=[]
    plotting_length=[]
    for record in SeqIO.parse(file_path, "fasta"):
        dna_sequence = str(record.seq).upper()
        special_count = 0
        total_leucine_serine_count = 0
        size_of_gene=len(dna_sequence)/3
        locus_tag = record.description.split()[0]
        upstream, downstream, beginner = find_flanking_sequences('GCA_000025985.1_ASM2598v1_genomic.fna.gbff', locus_tag)
        gc_content=calculate_gc_content(upstream + downstream)
        name=remove_first_word_preserve_whitespace(record.description)
        print(name)
        print(gc_content)
        print(size_of_gene)
        for i in range(0, len(dna_sequence), 3):
            codon = dna_sequence[i:i + 3]
            if codon in leucine_codons or codon in serine_codons:
                total_leucine_serine_count += 1
                if codon in special_codons:
                    special_count += 1

        if total_leucine_serine_count>=0:
            #proportion = (special_count - total_leucine_serine_count * 0.144)/np.sqrt(total_leucine_serine_count* 0.144*(1-0.144))
            #prob=stats.binom.cdf(special_count-1, total_leucine_serine_count,0.144)
            #proportion=np.log10(prob/(1-prob))
            #proportion=special_count/total_leucine_serine_count
            #gene_proportions[record.description+" "+str(gc_content)] = proportion
            #if "SusC" in record.description or "SusD" in record.description:
            #    susC_gene_proportions[record.description] = proportion
            plotting_deviations.append(special_count)
            to_translate=find_string_in_excel('out_mapped.xlsx',locus_tag)
            if len(to_translate)<5:
                to_translate=map_characters_to_entries(to_translate,cog_categories)
            else:
                print(to_translate)
                to_translate='Unknown'
            plotting_location.append(to_translate)
            print(to_translate)
            plotting_loci.append(find_string_in_excel('psortb_gramneg.xlsx',locus_tag+' '))
            plotting_gc.append(gc_content)
            plotting_name.append(name)
            plotting_start.append(beginner)
            plotting_size.append(total_leucine_serine_count)
            plotting_length.append(size_of_gene)
        if process_numeric_string(record.description) >= 4300:
            break
    # Determine the range and bins for the histograms
    # Extracting the values from the dictionaries
    all_genes_values = list(gene_proportions.values())
    #susC_genes_values = list(susC_gene_proportions.values())


    # Create a figure and a single subplot
    fig, ax1 = plt.subplots()
    #plotting_deviations_smooth=smooth_data_with_level(plotting_deviations,1)
    #plotting_gc_smooth=smooth_data_with_level(plotting_gc,1)
    clone_loci=[]
    clone_deviations=[]
    clone_gene=[]
    clone_gc=[]
    clone_size=[]
    for i in range(len(plotting_loci)):
        if True:#plotting_gc_smooth[i]<1 and plotting_gc_smooth[i]>0 :
            clone_loci.append(plotting_loci[i])
            clone_deviations.append(plotting_deviations[i])
            clone_gc.append(plotting_gc[i])
            clone_gene.append(plotting_name[i])
            clone_size.append(plotting_size[i])
    # Organize data into a dictionary
    data = {'Gene Name': plotting_name, 'COG Category': plotting_location, 'Cellular Location': plotting_loci, 'TTA/TCA Count': plotting_deviations,'Serine/Leucine Count':plotting_size,'Genome Position': plotting_start,'Flanking GC Content':plotting_gc,'Protein Length':plotting_length}

    # Create a DataFrame
    df = pd.DataFrame(data)

    # Output the DataFrame to an Excel file
    df.to_excel('genes_with_COG.xlsx', index=False, engine='openpyxl')


def calculate_gc_content_excluding_codons(dna_sequence):
    """
    Calculate the GC content of a DNA sequence (in lowercase), excluding the codons for leucine and serine.

    Args:
    dna_sequence (str): A string representing the DNA sequence.

    Returns:
    float: The GC content percentage, excluding leucine and serine codons.
    """
    # Convert to lowercase and define leucine and serine codons
    leucine_codons = {'tta', 'ttg', 'ctt', 'ctc', 'cta', 'ctg'}
    serine_codons = {'tca', 'tct', 'tcg', 'tcc', 'agc', 'agt'}

    excluded_codons = leucine_codons.union(serine_codons)

    # Initialize counters
    gc_count = 0
    total_count = 0

    # Iterate over the sequence in steps of 3 (codon size)
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i+3]
        if codon not in excluded_codons:
            gc_count += codon.count('g') + codon.count('c')
            total_count += 3  # Counting each nucleotide in the codon

    # Calculate GC content percentage
    if total_count == 0:
        return 0  # Avoid division by zero
    else:
        return (gc_count / total_count)


def process_numeric_string(s):
    # Find the index of the substring "BJOKPJ_"
    index = s.find("BJOKPJ_")
    if index == -1:
        return "The pattern BJOKPJ_ was not found in the input string."

    # Move the index to start after "BJOKPJ_"
    start_index = index + len("BJOKPJ_")

    # Extract the numeric characters immediately following "BJOKPJ_"
    numeric_part = ''
    for char in s[start_index:]:
        if char.isdigit():
            numeric_part += char
        else:
            # Stop at the first non-numeric character
            break

    # Attempt to convert the extracted numeric part to an integer
    if numeric_part:
        try:
            number = int(numeric_part)
        except ValueError:
            return "The numeric part following BJOKPJ_ could not be converted to an integer."

        # Divide the integer by 5
        result = number / 5

        return result
    else:
        return "No numeric part found following BJOKPJ_."


def find_flanking_sequences(genbank_file, locus_tag, flank_length=1000):
    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type == "gene" and "locus_tag" in feature.qualifiers:
                if feature.qualifiers["locus_tag"][0] == locus_tag:
                    feature_start = feature.location.start.position
                    feature_end = feature.location.end.position
                    strand = feature.location.strand

                    # Calculate upstream region
                    upstream_start = max(0, feature_start - flank_length)
                    upstream_end = feature_start
                    upstream_sequence = record.seq[upstream_start:upstream_end]
                    if strand == -1:  # Reverse strand
                        upstream_sequence = upstream_sequence.reverse_complement()

                    # Calculate downstream region
                    downstream_start = feature_end
                    downstream_end = min(len(record.seq), feature_end + flank_length)
                    downstream_sequence = record.seq[downstream_start:downstream_end]
                    if strand == -1:  # Reverse strand
                        downstream_sequence = downstream_sequence.reverse_complement()

                    return upstream_sequence, downstream_sequence, feature_start
    return None, None, None


file_path = "GCA_000025985.1_ASM2598v1_genomic.fasta"

run_main_code(file_path)
