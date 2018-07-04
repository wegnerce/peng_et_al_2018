'''
@author:      Carl-Eric Wegner
@affiliation: Küsel Lab - Aquatic Geomicrobiology
              Friedrich Schiller University of Jena

              carl-eric.wegner@uni-jena.de

Given a defined taxonomic group/level/, the below script extracts all the
sequences belonging to this group using multiple output files from QIIME1.9.1

** taxonomy assignments from assign_taxonomy.py [QIIME]
** picked OTUs [USEARCH]
'''

# modules to be imported
import csv
import time
from Bio import SeqIO

# files to be processed
# ** assigned_taxonomy, usearch_picked_otus, and sequence_data refer to output
#    from QIIME1.9.1 and usearch7, respectively
#    --> paths need to be modified accordingly
# ** the other variables are output files that are necessary/useful for downstream
#    analysis; collected_OTU_IDs (list of OTUs affiliated with the taxonomic group
#    of interest), collected_sequence_IDs (list of sequence IDs affiliated with the
#    taxonomic group of interest), collected_sequences (subsampled sequence data,
#    filtered for sequences affiliated with the taxonomic group of interest)

assigned_taxonomy = "/path/to/taxonomic.assignments"
collected_OTU_IDs = "/path/to/output/file/containing/OTU.ids"
usearch_picked_otus = "/path/to/usearch/cluster/file.uc"
collected_sequence_IDs = "/path/to/output/file/containing/OTU.ids"
sequence_data ="/output/from QIIMEs split_libraries.py script/seqs.fna"
collected_sequences = open("/path/output/file/seqs/of/interest.fna", "wb")

to_look_for = ["group_of_interest"] # this string needs to be modified dependent on
# the group you are interested in
IDs_of_interest = []
sequence_IDs_of_interest = []
sequences_of_interest = []

# STEP -I- check taxonomy file for OTU IDs of interest
with open(assigned_taxonomy) as infile:
    for row in csv.reader(infile, delimiter='\t'):
        for wanted in to_look_for:
            if wanted in row[1]:
                IDs_of_interest.append(row[0])
    infile.close()

with open(collected_OTU_IDs, 'wb') as outfile:
    collected_otus = csv.writer(outfile)
    for item in IDs_of_interest:
        print 'Collected ID: '+ item
        collected_otus.writerow ([item])
    outfile.close()

print '%i OTUs of interest detected and saved!' % len(IDs_of_interest)
time.sleep(10)

# STEP -II- extract sequence IDs from USEARCH picked_otus.py output file
with open(usearch_picked_otus) as infile:
    for row in csv.reader(infile, delimiter='\t'):
        if row[9] in IDs_of_interest:
            sequence_IDs_of_interest.append(row[8][:row[8].find(" ")])
    infile.close()

print '%i Sequence IDs of interest identified!' % len(sequence_IDs_of_interest)
time.sleep(10)

with open(collected_sequence_IDs, 'wb') as outfile:
    outfile.writelines("%s\n" % str(item) for item in sequence_IDs_of_interest)
    outfile.close()

print 'Sequence IDs saved... Going to extract sequences now!'
time.sleep(10)

# STEP -III- pick sequences of interest from the set of original sequences, meaning the output
# of QIIMEs split_libraries.py script
#wanted = set(line.rstrip("\n").split(None,1)[0] for line in open(collected_sequence_IDs))
print "Extract now %i sequences from the original data set." % (len(sequence_IDs_of_interest))

index = SeqIO.index(sequence_data, "fasta")
records = (index[r] for r in sequence_IDs_of_interest)
count = SeqIO.write(records, collected_sequences, "fasta")
assert count == len(sequence_IDs_of_interest)

print "Extracted %i sequences from %s (original data) to %s (taxon-specific subset)." % (count, sequence_data, collected_sequences)
