# -*- coding: utf-8 -*-
"""
@author:      Carl-Eric Wegner
@affiliation: Küsel Lab - Aquatic Geomicrobiology
              Friedrich Schiller University of Jena

              carl-eric.wegner@uni-jena.de

"""

###############################################################################
# NEEDED MODULES, PATHS, PARAMETERS
###############################################################################
import argparse
import sqlite3
import subprocess
import textwrap
from Bio import SeqIO

db = "dbcan.db"
no_annotated = 0
no_seqs = 0

seq_filter_ids = {"Cellulose utilization" : [], "Chitin utilization" : [],
                  "Xylan utilization" : [], "Pectin utilization" : [],
                  "Other CAZymes" : []}

func_filter = {}
func_seqs = {}

###############################################################################
# DEFINED CAZY MODULES
###############################################################################
cellulose_mapping = {"3.2.1.4_Cellulase" : 0, "3.2.1.91_Cellulose 1,4-beta-cellobiosidase" : 0,
                     "3.2.1.150_Oligoxyloglucan reducing-end-specific cellobiohydrolase" : 0,
                     "2.4.1.20_Cellobiose phosphorylase" : 0,
                     "1.1.99.19_Cellobiose dehydrogenase" : 0}

chitin_mapping = {"3.2.1.14_Chitinase" : 0, "3.5.1.41_Chitin deacetylase" : 0,
                  "3.5.1.21_N-acetylglucosamine-6-phosphate deacetylase" : 0,
                  "3.5.1.28_N-acetylmuramoyl-L-alanine amidase" : 0,
                  "3.5.1.89_N-acetylglucosaminylphosphatidylinositol deacetylase" : 0,
                  "3.2.1.52_N-acetylglucosaminidase" : 0}

xylan_mapping = {"3.2.1.8_Endo-1,4-beta-xylanase" : 0, "3.2.1.37_Xylan 1,4-beta xylosidase" : 0,
                 "3.1.1.72_Acetylxylan esterase" : 0,
                 "3.2.1.151_Xyloglucanase" : 0}

pectin_mapping = {"4.2.2.10_Pectin lyase" : 0, "3.1.1.11_Pectinesterase" : 0}

other_mapping = {"3.2.1.113_Mannosidase" : 0, "3.2.1.22_Alpha-galactosidase" : 0,
                 "3.2.1.23_Beta-galactosidase" : 0, "3.2.1.51_Alpha-fucosidase" : 0}

cellulose_ids = ("3.2.1.4_Cellulase", "3.2.1.91_Cellulose 1,4-beta-cellobiosidase",
                 "3.2.1.150_Oligoxyloglucan reducing-end-specific cellobiohydrolase",
                 "2.4.1.20_Cellobiose phosphorylase",
                 "1.1.99.19_Cellobiose dehydrogenase")

chitin_ids = ("3.2.1.14_Chitinase", "3.5.1.41_Chitin deacetylase",
              "3.5.1.21_N-acetylglucosamine-6-phosphate deacetylase",
              "3.5.1.28_N-acetylmuramoyl-L-alanine amidase",
              "3.5.1.89_N-acetylglucosaminylphosphatidylinositol deacetylase",
              "3.2.1.52_N-acetylglucosaminidase")

xylan_ids = ("3.2.1.8_Endo-1,4-beta-xylanase", "3.2.1.37_Xylan 1,4-beta xylosidase",
             "3.1.1.72_Acetylxylan esterase",
             "3.2.1.151_Xyloglucanase")

pectin_ids = ("4.2.2.10_Pectin lyase", "3.1.1.11_Pectinesterase")

other_ids = ("3.2.1.113_Mannosidase", "3.2.1.22_Alpha-galactosidase",
                 "3.2.1.23_Beta-galactosidase", "3.2.1.51_Alpha-fucosidase")

###############################################################################
# FUNCTIONS
###############################################################################

# parse Gbk-IDs for EC numbers
def check_dbCAN(gbk_id):
    conn = sqlite3.connect(db) # connect to our database
    cursor = conn.cursor()
    command = "SELECT ec FROM dbCAN WHERE family = '" + gbk_id +  "';" # lookup GO tags
    cursor.execute(command)
    results = list(set([item[0] for item in cursor.fetchall()]))
    cursor.close()
    if results: # results are retrieved as list
        #return result[0].encode('utf-8')
        return results
    else:
        return "no hit"

def seq_filter_function(diamond_out, functions):
    with open(functions, "rb") as infile:
        for line in infile:
            func_filter[line[:-2]] = []
    with open(diamond_out, "rb") as infile:
        for line in infile:
            look_up = check_dbCAN(line.split("\t")[1].split("|")[1])
            for entry in look_up: # the magic happens here
                if look_up != "no hit":
                    hit_key = next((k for (k,v) in func_filter.iteritems() if entry == k.split("_")[0]), None)
                    if hit_key:
                        for key in func_filter:
                            if hit_key in key:
                                func_filter[key].append(line.split("\t")[0])

def seq_filter(diamond_out):
    with open (diamond_out, "rb") as infile:
        # we look up CAZy families, however this time we want to identify sequences
        # that are linked to our defined CAZy modules
        for line in infile:
            look_up = check_dbCAN(line.split("\t")[1].split("|")[1])
            for entry in look_up: # the magic happens here
                if look_up != "no hit":
                    hit_key = next((k for (k,v) in cellulose_mapping.iteritems() if entry == k.split("_")[0]), None)
                    if hit_key:
                        seq_filter_ids["Cellulose utilization"].append(line.split("\t")[0])
                    hit_key = next((k for (k,v) in chitin_mapping.iteritems() if entry == k.split("_")[0]), None)
                    if hit_key:
                        seq_filter_ids["Chitin utilization"].append(line.split("\t")[0])
                    hit_key = next((k for (k,v) in xylan_mapping.iteritems() if entry == k.split("_")[0]), None)
                    if hit_key:
                        seq_filter_ids["Xylan utilization"].append(line.split("\t")[0])
                    hit_key = next((k for (k,v) in pectin_mapping.iteritems() if entry == k.split("_")[0]), None)
                    if hit_key:
                        seq_filter_ids["Pectin utilization"].append(line.split("\t")[0])
                    hit_key = next((k for (k,v) in other_mapping.iteritems() if entry == k.split("_")[0]), None)
                    if hit_key:
                        seq_filter_ids["Other CAZymes"].append(line.split("\t")[0])

# determine seq counts for defined functional modules
def annot(diamond_out, seqs):
    global no_annotated, no_seqs
    with open (diamond_out, "rb") as infile:
        no_annotated = subprocess.check_output("wc -l %s"%diamond_out, shell=True)
        no_annotated = no_annotated.split(" ")[0]
        no_seqs = len(list(SeqIO.parse(seqs, "fasta")))
        for line in infile:
            look_up = check_dbCAN(line.split("\t")[1].split("|")[1]) # we look up GO tags for Uniref
            if look_up != "no hit":
                for hit in look_up:
                    for k in cellulose_mapping:
                        if hit == k.split("_")[0]:
                            cellulose_mapping[k] = cellulose_mapping[k] + 1
                    for k in chitin_mapping:
                        if hit == k.split("_")[0]:
                            chitin_mapping[k] = chitin_mapping[k] + 1
                    for k in xylan_mapping:
                        if hit == k.split("_")[0]:
                            xylan_mapping[k] = xylan_mapping[k] + 1
                    for k in pectin_mapping:
                        if hit == k.split("_")[0]:
                            pectin_mapping[k] = pectin_mapping[k] + 1
                    for k in other_mapping:
                        if hit == k.split("_")[0]:
                            other_mapping[k] = other_mapping[k] + 1

# write output for respective program modes
def write_output(mode, diamond_out = None, seqs = None):
    if mode == "annot":
        with open(diamond_out + ".annot", 'w') as f:
            f.write("Cellulose utilization" + "\n")
            [f.write("{0},{1}\n".format(key, cellulose_mapping[key])) for key in cellulose_ids]
            f.write("Chitin utilization" + "\n")
            [f.write("{0},{1}\n".format(key, chitin_mapping[key])) for key in chitin_ids]
            f.write("Xylan utilization" "\n")
            [f.write("{0},{1}\n".format(key, xylan_mapping[key])) for key in xylan_ids]
            f.write("Pectin utilization" + "\n")
            [f.write("{0},{1}\n".format(key, pectin_mapping[key])) for key in pectin_ids]
            f.write("Other CAZymes" + "\n")
            [f.write("{0},{1}\n".format(key, other_mapping[key])) for key in other_ids]
            f.write("{0},{1}\n".format("No. of sequences with hits in CAZy", str(no_annotated)))
            f.write("{0},{1}\n".format("No. of sequences", str(no_seqs)))
    elif mode == "filter":
        cellulose_seqs = (r for r in SeqIO.parse(seqs, "fasta") if r.id in seq_filter_ids["Cellulose utilization"])
        count = SeqIO.write(cellulose_seqs, seqs + ".cellulose.fa", "fasta")
        print " *** Saved %i sequences linked to Cellulose utilization." % (count)
        chitin_seqs = (r for r in SeqIO.parse(seqs, "fasta") if r.id in seq_filter_ids["Chitin utilization"])
        count = SeqIO.write(chitin_seqs, seqs + ".chitin.fa", "fasta")
        print " *** Saved %i sequences linked to Chitin utilization." % (count)
        xylan_seqs = (r for r in SeqIO.parse(seqs, "fasta") if r.id in seq_filter_ids["Xylan utilization"])
        count = SeqIO.write(xylan_seqs, seqs + ".xylan.fa", "fasta")
        print " *** Saved %i sequences linked to Xylan utilization." % (count)
        pectin_seqs = (r for r in SeqIO.parse(seqs, "fasta") if r.id in seq_filter_ids["Pectin utilization"])
        count = SeqIO.write(pectin_seqs, seqs + ".pectin.fa", "fasta")
        print " *** Saved %i sequences linked to Pectin utilization." % (count)
        other_seqs = (r for r in SeqIO.parse(seqs, "fasta") if r.id in seq_filter_ids["Other CAZymes"])
        count = SeqIO.write(other_seqs, seqs + ".other.fa", "fasta")
        print " *** Saved %i sequences linked to Other CAZymes." % (count)
    elif mode == "func_filter":
        for k, v in func_filter.iteritems():
            func_seqs = (r for r in SeqIO.parse(seqs, "fasta") if r.id in func_filter[k])
            count = SeqIO.write(func_seqs, seqs + "." + k + ".fasta", "fasta")
            print" Saved %i sequences linked to %s." % (count, k)

###############################################################################
# MAIN PROGRAM STRUCTURE
###############################################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog = "dbCAN_annot.py", # argumment parsing
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = textwrap.dedent ('''\
    | dbCAN ANNOTATOR                                        |
    | v0.1 November, '17                                     |
    | ------------------------------------------------------ |
    |                                  (c) Carl-Eric Wegner  |
    |                       for comments, reporting issues:  |
    |                          carl-eric.wegner@uni-jena.de  |
    '''))
    parser.add_argument("--mode", "-m", help="annotate ('annot') based on diamond output or filter ('filter') sequences module-based", action="store")
    parser.add_argument("--diamond", "-d", help="diamond output to be processed (.tab formatted)", action="store")
    parser.add_argument("--seq", "-s", help="sequence data queried against dbCAN (.fasta)", action="store" )
    parser.add_argument("--func", "-f", help="enzyme functions of interest for sequence filtering", action="store" )
    args = parser.parse_args() # we need only diamond tab separated output as essential input
    print ""
    print " | dbCAN ANNOTATOR                                        | "
    print " | v0.1 November, '17                                     |."
    print " | ------------------------------------------------------ | "
    print " |                                  (c) Carl-Eric Wegner  | "
    print " |                       for comments, reporting issues:  | "
    print " |                          carl-eric.wegner@uni-jena.de  | "
    print ""
    if args.diamond and args.seq and args.mode == "annot":
        print " *** Determine Annotation statistics for defined CAZY modules."
        annot(args.diamond, args.seq)
        print " *** Write annotation output."
        write_output(args.mode, args.diamond)
        print ""
        print " *** Done."
    elif args.diamond and args.seq and args.mode == "filter":
        print " *** Filter sequences based on defined CAZy modules."
        seq_filter(args.diamond)
        print " *** Write filtered sequences."
        write_output(args.mode,None,args.seq)
        print ""
        print " *** Done."
    elif args.diamond and args.seq and args.func and args.mode == "func_filter":
        print " *** Filter sequences based on provided list of enzyme functions."
        seq_filter_function(args.diamond, args.func)
        print " *** Write filtered sequences."
        write_output(args.mode,None,args.seq)
        print ""
        print " *** Done."
    else:
        print " *** Please specify necessary input for annotation/filtering:"
        print "       - output from diamond/(ublast) queries against dbCAN/CAZy "
        print "         tab-separated standard blast output (outfmt 6) "
        print "       - sequence data that was queried against dbCAN/CAZy"
        print "         in .fasta format"
        print "         --> mode 'annot' for annotation, 'filter' for filtering"
