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
    '-d', dest='d', help='Vibrio,all', required=True)
parser.add_argument(
    '-n', dest='n', help='selected genus', default="Vibrio")
parser.add_argument(
    '-x', dest='x', help='excluded genus', default="")


args = parser.parse_args()

def checking_copy1(linea):
    copia=False
    # From genus Vibrio: Total number of strains 100% identifiable:
    # 109 strains from 47 species
    hits = linea.split(":")[2]
    sts = int(hits.split("strains")[0].replace(" ", ""))
    sps = hits.split("from")[1].replace("species", "")
    sps = int(sps.replace(" ", ""))
    if sts > int(args.t) and sps > int(args.s):
        copia = True
    return copia

def checking_copy2(lin):
    copia2 = True
    copia = True
    list_exclu=[l.capitalize() for l in args.x.split(",") if l != ""]
    if len(list_exclu) > 0 :
        exclu=lin.split(":")[1]
        exclu=[e.lstrip() for e in exclu.split(",")]
        for item in exclu:
            if item in list_exclu:
                copia2 = False
                copia = False
                break
    return copia, copia2

def print_out(h,l,f1,f2, tp):
    if "---" in l:
        text= "*** full Genus specificity"
    else:
        text= l[1:]
    print(h, "\n", text, file=f1)
    print(" ".join(tp), file=f2)
    copia=False
    copia2=False
    return copia, copia2

file_list = args.o + "_list.txt"

with open(args.i, "r") as fin, open(args.o, "w") as fout, open(file_list, "w") as flist:
    print("#Gene_target Forward_primer sequence_(5'-> 3') Reverse_primer sequence_(5'-> 3')", file=flist)
    copy2=False
    copy=False
    for line in fin:
        line = line.rstrip()
        if line.startswith("#Target"):
            # Target gene: GroEL_124 - forward primer: Primer_18_at_95
            # YCCNAAAGGHCGYAAYG - reverse primer: Primer_178_at_567
            # AAYTGCATVCCYTCNAC
            header = line
            toprint = [line.split()[i] for i in [2, 6, 7, 11, 12]]

        if args.d.capitalize() == args.n.capitalize():
            str="From genus "+args.n.capitalize()+": Total number of strains"
            if str in line:
                copy=checking_copy1(line)
                if copy:
                    print(header, "\n", line[1:], file=fout)
                    print(" ".join(toprint), file=flist)
                    copy=False

        else:

            str2="Number of amplicons that are not "+args.n.capitalize()
            strin="where, genus (exluding"
            str3=" -------------------------------------"

            if str2 in line:
                    # Number of amplicons that are not Vibrio: Genus 303 Spp
                    # 690 strains 691
                hits = line.split(":")[1]
                gns = hits.split("Spp")[0].replace("Genus", "")
                gns = int(gns.replace(" ", ""))
                if gns <= int(args.g):
                    copy = True
            if copy and (str3 in line):
                copy2=True
            if copy and (strin in line):
                copy, copy2=checking_copy2(line)
            if copy and copy2:
                copy, copy2=print_out(header, line, fout, flist, toprint)
