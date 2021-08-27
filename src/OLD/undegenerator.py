#!/usr/bin/env python3

import os
import argparse
#import NumPy as np


usage = 'undegenerator.py -f fwd -F fwdname -r rev -R revname'
description = 'This program creates a list of primers from degenerated primers'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument(
    '-f',
    dest='pf',
    help='degenerated forward primer sequence')
parser.add_argument(
    '-r',
    dest='pr',
    help='degenerated reverserse primer sequence')
parser.add_argument('-F', dest='outf', help='Forward primer name')
parser.add_argument('-R', dest='outr', help='Reverse primer name')
parser.add_argument('-o', dest='o', help='Output directory')

args = parser.parse_args()

if args.outf and args.outr:
    both = args.outf + "_" + args.outr

degen = {"R": ["A", "G"], "Y": ["C", "T"], "M": ["A", "C"],
         "S": ["G", "C"], "W": ["A", "T"], "K": ["G", "T"],
         "V": ["A", "G", "C"], "D": ["A", "G", "T"],
         "H": ["A", "T", "C"], "B": ["C", "G", "T"],
         #       "I": "N" }
         "I": ["A", "G", "C", "T"], "N": ["A", "G", "C", "T"]}

#RYMKSWVBDHIN / YRKMSWBVHDNN
if not os.path.exists(args.o):
    os.makedirs(args.o)


def undegenerating(primer, degen, namep, type):
    total_deg = 1
    changes = []  # List of changes per position
    for item in primer:
        if item in degen.keys():
            total_deg = total_deg * len(degen[item])  # Total number of primers
            changes.append(total_deg)
        else:
            changes.append(1)

    primers = {}  # Dictionary storing the primers

    for i in range(0, total_deg):
        primers[i] = ["N" for n in range(0, len(primer))]  # initializing

    for it in range(0, len(primer)):
        deg = primer[it]
        pr = 0
        while pr < total_deg:  # for all the primers
            if deg in degen.keys():
                Nn = len(degen[deg])  # number of nucleotides
                Ntc = changes[it]  # number of changes
                # number of times the same nucleotide should be used for if
                # window change
                window = total_deg / Ntc
                for ch in range(
                        0, Nn):  # for each nucleotide in degenerated nucleotide
                    counter = 0
                    while counter < window:
                        primers[pr][it] = degen[deg][ch]
                        counter += 1
                        pr += 1
            else:
                primers[pr][it] = deg
                pr += 1

    # unique_list=[]
    with open(os.path.join(args.o, namep + type + ".fasta"), "w") as ffasta, open(os.path.join(args.o, both + ".fasta"), "a") as fboth:
        for k in primers.keys():
            #    print("".join(primers[k]), file=fout)
            print(">primer_{}_{}_{}\n{}".format(k, namep, type, "".join(primers[k])), file=ffasta)
            print(">primer_{}_{}_{}\n{}".format(k, namep, type, "".join(primers[k])), file=fboth)
    print(
        "        total number of unique sequences from {} {}: {}".format(
            namep,
            type,
            total_deg))


    #        unique_list.append("".join(primers[k]))
# print(len(unique_list))
# print(len(set(unique_list)))
if (args.pf and args.outf):
    undegenerating(args.pf, degen, args.outf, "forward")
else:
    print("Forward primer information no provided")

if (args.pr and args.outr):
    undegenerating(args.pr, degen, args.outr, "reverse")
else:
    print("Reverse primer information no provided")
