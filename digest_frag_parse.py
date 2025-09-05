"""
digest_frag_parse.py
Created on June 6th, 2019 by Myung Chang Lee (Noah Lee)
Parses text result from https://www.bioinformatics.org/sms2/rest_digest.html
Restriction digest output

"""

# Import sys to get arguments
import sys
# Import re for regex
import re

# Get the first argument from command line, filename
try:
    in_filename = sys.argv[1]
except Exception as e:
    # If no filename was supplied, then use the default name
    in_filename = input("Please enter the filename: ")

if len(sys.argv[1:]) != 1:
    print("Please do note that you can specify input and filename on commandline")
    print("E.g. python3 digest_frag_parse.py digest_result.txt")

with open(in_filename, "rt") as in_file:
    fragments = dict()

    for line in in_file.readlines():
        match_obj = re.search(r">(\d+) bp.+parent (.+),.+", line)
        if match_obj:
            if match_obj.group(2) in fragments:
                fragments[match_obj.group(2)].append(match_obj.group(1))
            else:
                fragments.update({match_obj.group(2): [match_obj.group(1)]})

    # sort plasmids alphabetically
    for plasmid in sorted(fragments, key=str.lower):
        print("For", plasmid, ":")
        print(" ----  ", "bp, ".join(fragments[plasmid]) + "bp")
