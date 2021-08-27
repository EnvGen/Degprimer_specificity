#!/usr/bin/env python3

import os
import argparse

usage = 'python parse_Summary.py -i -o -s -t -d -g -n'
description = 'This program selects primers according to amplicon characteristics'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument(
    '-i', dest='i', help='input file Summary.txt', required=True)
parser.add_argument(
    '-o', dest='o', help='output file selected_primer.txt', required=True)
parser.add_argument(
    '-s', dest='s', help='minimum number of identifiable species', default=48)
parser.add_argument(
    '-t', dest='t', help='minimum number of identifiable strains', default=153)
parser.add_argument(
    '-g', dest='g', help='maximum number of identifiable Genus', default=50)
parser.add_argument(
    '-d', dest='d', help='Vibrio,Bacteria', required=True)
parser.add_argument(
    '-n', dest='n', help='selected genus', default="Vibrio")

args = parser.parse_args()


file_list = args.o + "_list.txt"


with open(args.i, "r") as fin, open(args.o, "w") as fout, open(file_list, "w") as flist:
    print("#Gene_target Forward_primer sequence_(5'-> 3') Reverse_primer sequence_(5'-> 3')", file=flist)
    for line in fin:
        line = line.rstrip()
        if line.startswith("#Target"):
            # Target gene: GroEL_124 - forward primer: Primer_18_at_95
            # YCCNAAAGGHCGYAAYG - reverse primer: Primer_178_at_567
            # AAYTGCATVCCYTCNAC
            header = line
            toprint = [line.split()[i] for i in [2, 6, 7, 11, 12]]

        if args.d == args.n:
            str="From genus "+args.n+": Total number of strains"
            if str in line:
                # From genus Vibrio: Total number of strains 100% identifiable:
                # 109 strains from 47 species
                hits = line.split(":")[2]
                sts = int(hits.split("strains")[0].replace(" ", ""))
                sps = hits.split("from")[1].replace("species", "")
                sps = int(sps.replace(" ", ""))
                if sts > int(args.t) and sps > int(args.s):
                    print(header, "\n", line[1:], file=fout)
                    print(" ".join(toprint), file=flist)

        else:
            str2="Number of amplicons that are not "+args.n
            if str2 in line:
                    # Number of amplicons that are not Vibrio: Genus 303 Spp
                    # 690 strains 691
                hits = line.split(":")[1]
                gns = hits.split("Spp")[0].replace("Genus", "")
                gns = int(gns.replace(" ", ""))
                if gns <= int(args.g):
                    print(header, "\n", line[1:], file=fout)
                    print(" ".join(toprint), file=flist)
