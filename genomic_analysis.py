# The first thing to do is open our file with the downloaded FASTA data file.
# https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/
data = open('Homo_sapiens.GRCh38.dna.chromosome.21.fa', 'r')
# Read the contents of the file
dna_seq = data.read()
print("DATA FILE: \n" + dna_seq[:1500] + "\n")
print("DNA SEQUENCE LENGTH:" + str(len(dna_seq)) + "\n")
# Extract the FASTA header to separate it out from the genetic code.
fasta_header = dna_seq.partition('\n')[0]
print("FASTA HEADER" + fasta_header)

# Remove the FASTA header from the sequence data.
print("HEADER LENGTH = " + str(len(fasta_header)) + "\n")
dna_seq = dna_seq[56:]

seq_len = len(dna_seq) # DNA sequence length
print("FASTA DATA SEQUENCE LENGTH = " + str(seq_len) + "\n")

# Find positions in the middle of the sequence which contain 'N'.
for pos in range (1, len(dna_seq) - 1):
  if dna_seq[pos] == 'N':
    if dna_seq[pos - 1] != 'N' and dna_seq[pos + 1] != 'N':
      print(pos)

# Function that finds the start of the well-identified coding region (A, C, G, 
  # and T instead of N)
def find_start_coding_region(seq):
  for pos in range(1, len(seq) - 1):
    # We need to check if the position two bases out from the current position
      # is also not 'N' because there could be a newline character in the following
      # position.
    if seq[pos] == 'N' and seq[pos + 1] != 'N' and seq[pos + 2] != 'N':
      return pos + 1

# Function that find the position of the last chunk of unidentified region in 
  # the sequence (long string of N's)
def find_end_coding_region(seq):
  last_N_chunk = 0
  for pos in range(1, len(seq) - 1):
    # Similarly, we need to check whether two positions before the current position
      # is not an 'N'.
    if seq[pos] == 'N' and seq[pos - 1] != 'N' and seq[pos - 2] != 'N':
      last_N_chunk = pos
  return last_N_chunk - 1 # Return the last occurence (at the very end of the sequence)

# Store these start and end positions (the sequence in between we'll call our 
  # "coding region")
start_coding_region = find_start_coding_region(dna_seq)
print("Start coding region: " + str(start_coding_region))

end_coding_region = find_end_coding_region(dna_seq)
print("End coding region: " + str(end_coding_region))

print("Length of N's at end :" + str(seq_len - end_coding_region))

print("Sequence length remaning: " + str(seq_len - start_coding_region - (seq_len - end_coding_region)))

# Since the sequence contains over 47 million base pairs, trying to parse all
  # the data at once requires too much processing power and results in a runtime 
  # disconnect error. To work around this issue, we'll parse the sequence in 
  # chunks of four, of size 10 million base pairs in each iteration.

# Start and end position for sequence parse in first iteration of loop
start_pos = start_coding_region
end_pos = 10000000

# We'll build up a new sequence of the modified DNA, with all '\n removed.
mod_dna_seq = ""

# Modify DNA sequence to remove all \n' chars
print("Modifying DNA Sequence")
while end_pos < end_coding_region:
  dna_chunk = dna_seq[start_pos: end_pos].replace('\n', '') # Remove each instance of '\n'
  mod_dna_seq += dna_chunk # Build up modified DNA sequence
  # Print the start and end positions in each iteration to make sure we're parsing
    # all segments of the sequence
  print("start: " + str(start_pos) + " end: " + str(end_pos) + " dna length: " +
        str(len(mod_dna_seq)))
  start_pos = end_pos + 1 # Update start position
  end_pos += 10000000 # Increment end position by 10,000,000

# We still need to modify the last chunk of the FASTA file, up to the end of the
  # sequence, and then we'll have our final modified DNA sequence.
last_dna_chunk = dna_seq[start_pos: end_coding_region].replace('\n', '') # Remove each instance of '\n'
mod_dna_seq += last_dna_chunk # Build up modified DNA sequence
print("start: " + str(start_pos) + " end: " + str(seq_len) + " dna length: " + 
      str(len(mod_dna_seq)))

dna_length = len(mod_dna_seq) # Modified DNA sequence length
print("\nDNA SEQUENCE MODIFIED LENGTH: " + str(dna_length) + "\n")

# Difference in sequence length between original FASTA data and modified sequence
print("ELIMINATED SEQUENCE: " + str(seq_len - dna_length) + "\n")

# TRANSCRIPTION (DNA to mRNA)
mrna_seq = mod_dna_seq.replace('T', 'U')
mrna_length = len(mrna_seq)

print("mRNA LENGTH = " + str(mrna_length) + "\n")
print("SEQUENCES TRANSCRIBED CORRECTLY: " + str(dna_length == mrna_length) + "\n")

mrna_seq[:1500]

# Let's begin by finding the first occurence of the start codon in the DNA.
start_codon = mrna_seq.find('AUG')
print("Start codon: " + str(start_codon))

# We'll create a dictionary to encode the codons corresponding to each of the 20
  # amino acids, the start codon, the stop codons, and deal with N's in the sequence
  # as well.
codons = {"GCG" : "A", "GCA" : "A", "GCC" : "A", "GCU" : "A", # Alanine
            "UGC" : "C", "UGU" : "C", # Cysteine
            "GAC" : "D", "GAU" : "D", # Aspartic acid
            "GAG" : "E", "GAA" : "E", # Glutamic acid
            "UUC" : "F", "UUU" : "F", #Phenylalanine
            "GGG" : "G", "GGA" : "G", "GGC" : "G", "GGU" : "G", # Glycine
            "CAU": "H", "CAC" : "H", # Histidine
            "AUU" : "I", "AUC" : " I", "AUA" : "I", # Isoleucine
            "AAA" : "K", "AAG" : "K", # Lysine
            "CUU" : "L", "CUC" : "L", "CUA" : "L", "CUG" : "L", "UUG" : "L", "UUA" : "L", # Leucine
            "AUG" : "M", # Methionine (start codon)
            "AAC" : "N", "AAU" : "N", # Asparagine
            "CCU" : "P", "CCC" : "P", "CCA" : "P", "CCG" : "P", # Proline
            "CAA" : "Q", "CAG" : "Q", # Glutamine
            "CGU" : "R", "CGC" : "R", "CGA" : "R", "CGG" : "R", "AGA" : "R", "AGG" : "R", # Arginine
            "UCG" : "S", "UCA" : "S", "UCC" : "S", "UCU" : "S", "AGU" : "S", "AGC" : "S", # Serine
            "ACU" : "T", "ACC" : "T", "ACA" : "T", "ACG" : "T", # Threonine
            "GUG" : "V", "GUA" : "V", "GUC" : "V", "GUU" : "V", # Valine
            "UGG" : "W", # Tryptophan
            "UAC" : "Y", "UAU" : "Y", # Tyrosine
            "UGA" : "*", "UAG" : "*", "UAA" : "*", # Stop codons
          
            # Since we can't remove N's from the middle of the sequence (would result
              # in frameshift mutations), we'll enocde any codon in the mRNA that
              # has an N as "X" in the protein sequence. Here are all the possible
              # configurations in which N can be present in a codon:
            "NNN" : "X", 
            "NNA" : "X", "NNC" : "X", "NNG" : "X", "NNU" : "X",
            "ANN" : "X", "CNN" : "X", "GNN" : "X", "UNN" : "X",
            "NAA" : "X", "NAC" : "X", "NAG" : "X", "NAU" : "X",
            "NCA" : "X", "NCC" : "X", "NCG" : "X", "NCU" : "X",
            "NGA" : "X", "NGC" : "X", "NGG" : "X", "NAU" : "X",
            "NUA" : "X", "NUC" : "X", "NUG" : "X", "NUU" : "X",
            "ANA" : "X", "ANC" : "X", "ANG" : "X", "ANU" : "X",
            "CNA" : "X", "CNC" : "X", "CNG" : "X", "CNU" : "X",
            "GNA" : "X", "GNC" : "X", "GNG" : "X", "GNU" : "X",
            "UNA" : "X", "UNC" : "X", "UNG" : "X", "UNU" : "X",
            "AAN" : "X", "ACN" : "X", "AGN" : "X", "AUN" : "X",
            "CAN" : "X", "CCN" : "X", "CGN" : "X", "CUN" : "X",
            "GAN" : "X", "GCN" : "X", "GGN" : "X", "GUN" : "X",
            "UAN" : "X", "UCN" : "X", "UGN" : "X", "UUN" : "X",
            }

# TRANSLATION (mRNA to protein)
# This function synthesizes a specific protein sequence, starting from the start
  # codon in the mRNA sequence and terminating translation once a stop codon is reached.
def translate(start):
  protein_seq = "" # We'll build up the protein sequence
  translate_pos = start
  while mrna_seq[translate_pos: translate_pos+3] not in stop_codons: # Tranverse mRNA sequence by triplets
    amino_acid = codons[mrna_seq[translate_pos: translate_pos+3]] # Identify each amino acid based on codon dict
    protein_seq += amino_acid # Add the amino acid to the protein sequence
    translate_pos += 3
  protein_seq += amino_acid # Add stop codon to sequence
  proteins.append(protein_seq) # Add synthesized protein to proteins list
  return translate_pos # Return position of stop codon

# Now we'll call the translate function on our sequence
proteins = [] # Create a list to store all synthesized proteins
start_codon_seq = 'AUG' # Start codon
stop_codons = ['UGA', 'UAG', 'UAA'] # Stop codons

# Begin translation at the first start codon
translate_pos = translate(start_codon)

# Iterate through mRNA sequence, calling translate() each time start codon is found
while translate_pos < end_coding_region:
  if mrna_seq[translate_pos: translate_pos + 3] == start_codon_seq:
    translate_pos = translate(translate_pos)
  else:
    translate_pos += 3

print("Ended Translating at: " + str(translate_pos) + "\n")
print("PROTEINS: ")
proteins[:25] # Display first 25 proteins in list

import pandas as pd

# Let's start by creating a dataframe in which to store the protein seqeunces we've extracted
ch21code_df = pd.DataFrame(data=proteins, columns=['amino acid sequence']) # Store the amino acid sequences of each protein
ch21code_df['protein length'] = ch21code_df.apply(lambda row: len(row['amino acid sequence']), axis=1) # Store the amino acid lengths
ch21code_df['nucleotides'] = ch21code_df['protein length'] * 3 # Store the nucleotide lengths
print("Number of Proteins i DataFrame: " + str(len(ch21code_df)))
ch21code_df