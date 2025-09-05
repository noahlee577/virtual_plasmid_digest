"""
digest_frag.py
Created on June 6th, 2019 by Myung Chang Lee (Noah Lee)
Got tired of parsing digest results one by one so automated it
Input REs and get cut band sizes

"""

# Import sys to get arguments
import sys
# Import re for regex
import re

# Import restriction enzyme package
from Bio import Restriction

from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import seaborn as sns

# Get all the NEB enzymes
NEB_enzymes = Restriction.RestrictionBatch(suppliers=['N']).elements()

in_filename = "fasta.txt"

verbose = False
verbose_args = ["-v", "--verbose"]

# Get the first argument from command line, filename
try:
    in_filename = sys.argv[1]

    if len(sys.argv[1:]) <= 1:
        print("Please do note that you can specify filename and enzymes on commandline")
        print("E.g. python3 digest_frag_parse.py [" + "|".join(verbose_args) + "] fasta.txt enzymeI enzymeII...")

    if in_filename in verbose_args:
        verbose = True
        in_filename = sys.argv[2]
        enzymes = sys.argv[3:]
    else:
        enzymes = sys.argv[2:]

except Exception as e:
    # If no filename was supplied, then use the default name
    print("Please specify input filename and enzyme list")
    print("E.g. python3 digest_frag_parse.py fasta.txt EcoRI NcoI...")

for enzyme in enzymes:
    if enzyme in verbose_args:
        verbose = True
        enzymes.remove(enzyme)

    elif enzyme not in Restriction.AllEnzymes:
        print("Unrecognized enzyme:", enzyme)
        print("Please use proper capitalization/Roman numeral rules")
        print("ecori or EcoR1 != EcoRI")
        enzymes.remove(enzyme)

if not enzymes:
    enzymes = ["XmaI", "NcoI"]

amb = IUPACAmbiguousDNA()

plasmids = dict()

# The file needs to be formatted as the following:
# >Plasmid name1
# ATCG....whole sequence in one line
# >Plasmid name2
# ....
with open(in_filename, "rt") as in_file:
    # Put outside just so we can store them for later
    plasmid_name = "N/A"
    for line in in_file.readlines():
        match_obj = re.search(r">(.+)", line)
        if match_obj:
            plasmid_name = match_obj.group(1)
        else:
            plasmids[plasmid_name] = Seq(line, amb)

# Less efficient, but will create two separate ones
# all_data = Order by plasmid
# all_plasmid = Order by enzyme
all_data = dict()
all_plasmid = dict()

# For each plasmid, store restriction enzyme cut sequences produced by each enzyme
for plasmid in plasmids:
    RE_seqs = dict()
    RE_combo = dict()

    visited = list()

    # TODO: IMPLEMENT USING RECURSION

    for enzyme1 in ["XmaI"]:
        RE_seqs[enzyme1] = getattr(Restriction, enzyme1).catalyze(
            plasmids[plasmid], linear=False)
        if not getattr(Restriction, enzyme1).search(plasmids[plasmid], linear=False):
            # If not cut, the plasmid is still circular;
            # take this into account for secondary digestion
            print("No cutsite for ", enzyme1, "in plasmid", plasmid)
            circular = True
        else:
            circular = False

        for enzyme2 in enzymes:
            if enzyme1+enzyme2 in visited or enzyme2+enzyme1 in visited:
                continue

            visited.append(enzyme2+enzyme1)

            if enzyme1 == enzyme2:
                continue

            combo_key = enzyme1+"+"+enzyme2
            RE_combo[combo_key] = list()

            for first_cut_seq in RE_seqs[enzyme1]:
                if circular:
                    RE_combo[combo_key].append(
                        getattr(Restriction, enzyme2).catalyze(first_cut_seq, linear=False))
                else:
                    RE_combo[combo_key].append(
                        getattr(Restriction, enzyme2).catalyze(first_cut_seq, linear=True))

    # Get single digestion fragment lengths
    RE_seqs_len = dict()
    for key in RE_seqs:
        RE_seqs_len[key] = list()
        for element in RE_seqs[key]:
            RE_seqs_len[key].append(len(element))

    # Get double digest fragment lengths
    RE_combo_len = dict()
    for key in RE_combo:
        RE_combo_len[key] = list()
        for element in RE_combo[key]:
            for sequence in element:
                RE_combo_len[key].append(len(sequence))

    # Sort enzymes alphabetically
    # Also sort both fragments by descending fragment size
    RE_seqs_len = [(key, sorted(RE_seqs_len[key], reverse=True))
                   for key in sorted(RE_seqs_len)]
    RE_combo_len = [(key, sorted(RE_combo_len[key], reverse=True))
                    for key in sorted(RE_combo_len)]

    if verbose:
        print("For", plasmid + ":")
        print("---- ", RE_seqs_len)
        print("---- ", RE_combo_len)

    # for single_digest in RE_seqs_len:
    #     all_data[plasmid+" "+single_digest[0]] = single_digest[1]

    for double_digest in RE_combo_len:
        all_data[plasmid+" "+double_digest[0]] = double_digest[1]
        if double_digest[0] in all_plasmid:
            all_plasmid[double_digest[0]].update({plasmid: double_digest[1]})
        else:
            all_plasmid[double_digest[0]] = {plasmid: double_digest[1]}

plasmid_data_graphing = dict()

for enzyme_combo in all_plasmid:
    for plasmid, sizes in all_plasmid[enzyme_combo].items():
        plasmid_data_graphing[plasmid+" "+enzyme_combo] = sizes

if verbose:
    print(all_data)
    print(plasmid_data_graphing)

NEB_Ladder = [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500,
              1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100]
Ladder = {"NEB 1kb Plus": NEB_Ladder}

all_data.update(Ladder)
plasmid_data_graphing.update(Ladder)

# Convert to pandas dataframe
all_df = pd.DataFrame.from_dict(all_data, orient='index')
plasmid_df = pd.DataFrame.from_dict(plasmid_data_graphing, orient='index')

# Make all the "data points" a line to correspond to a band...
markers = ["_" for key in all_df]

if verbose:
    print(all_df)
    print(plasmid_df)

plasmid_df.T.to_excel("./output/Digest.xlsx")

#plt.figure(1, figsize=(8, 8))
#graph = sns.scatterplot(data=all_df, markers=markers, legend=False, s=200)
#sns.set(font_scale=0.5)
#graph.set(yscale="log")
#graph.figure.subplots_adjust(bottom=0.3)
#plt.xticks(rotation=90)

plt.figure()
graph2 = sns.scatterplot(data=plasmid_df, markers=markers, legend=False, s=200)
sns.set(font_scale=0.5)
#graph2.set(yscale="log")
graph2.figure.subplots_adjust(bottom=0.3)
plt.yscale('symlog', linthreshy=1000)
graph2.set_yticks(NEB_Ladder)
graph2.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.xticks(rotation=90)
plt.savefig("./output/20190605 DoubleDigests.png", dpi=600)
plt.show()
