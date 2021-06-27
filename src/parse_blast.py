#!/usr/bin/env python3

import os
import argparse
import re

usage = 'parse_blast.py -i -n -F -R -f -r -o -d -l'
description = 'This program counts the number of unique hits from a blast table, above an user-defined threshold'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument('-i', dest='inf', help='input file .blastn', required=True)
parser.add_argument('-n', dest='n', help='Target gene', required=True)
parser.add_argument('-F', dest='nf', help='Forward primer name', required=True)
parser.add_argument('-R', dest='nr', help='Reverse primer name', required=True)
parser.add_argument(
    '-f',
    dest='f',
    help='Original primer Forward sequence',
    required=True)
parser.add_argument(
    '-r',
    dest='r',
    help='Original primer reverse sequence',
    required=True)
parser.add_argument(
    '-d',
    dest='d',
    help='identity threshold default 99',
    default=99)
parser.add_argument(
    '-l',
    dest='l',
    help='percentage query coverage threshold default 100',
    default=100)
parser.add_argument('-o', dest='out', help='output file')
args = parser.parse_args()

def get_info(line):
    line = line.split()
   # id=line[1].split("_")[0] #it doesn't work since same strain may
   # have two different AccNumber because there are two Chromosomes
    id = line[1].split("_")
    id = "_".join(id[2:])
    id = re.sub(".chromosome.*$", "", id, count=1)
    id = re.sub(".genome.*$", "", id, count=1)
    id = re.sub(".plasmid.*$", "", id, count=1)
    id = re.sub(".DNA.*$", "", id, count=1)
    id = re.sub(".complete.*$", "", id, count=1)
    id = re.sub("sp.", "sp", id, count=1)
    bitscore = float(line[11])
    identity = float(line[2])
    align_len = float(line[3])
    primername = line[0]
    if re.search("forward", primername):
                primer_len = len(args.f)
    if re.search("reverse", primername):
                primer_len = len(args.r)
    p_align = round(align_len / primer_len * 100, 2)
    to_add = line[:3]
    to_add.append(str(primer_len))
    to_add.append(str(p_align))
    to_add.extend(line[3:])

    return id, bitscore, identity, to_add

def storing_blast_hits(file):
    id_hit = {}
    with open(file, "r") as fin:
        for line in fin:
            line = line.rstrip()
            if not line.startswith("#"):
                id, bitscore, identity, to_add=get_info(line)

                if id in id_hit:
                    if identity > float(id_hit[id][2]): #selecting the best hit
                        id_hit[id] = to_add
                    elif bitscore > float(id_hit[id][13]):
                        id_hit[id] = to_add
                else:
                    id_hit[id] = to_add
    return id_hit

def build_dict(key, mydict_to_add):
    g_s = key.split("_")[0]
    sp_s = key.split("_")[1]
    st_s = "_".join(key.split("_")[2:])
    if g_s in mydict_to_add:
        if sp_s in mydict_to_add[g_s].keys():
            if st_s in mydict_to_add[g_s][sp_s].keys():
                mydict_to_add[g_s][sp_s][st_s] += 1
            else:
                mydict_to_add[g_s][sp_s][st_s] = 1
        else:
            mydict_to_add[g_s][sp_s] = {st_s: 1}
    else:
        mydict_to_add[g_s] = {sp_s: {st_s: 1}}
    return mydict_to_add

def get_hits_by_genus(mydict, fout, fout2):
    my_selected = 0
    Genus = {}
    Genus_slc = {}
    for k in mydict.keys():
        Genus=build_dict(k, Genus)

        if (float(mydict[k][2]) >= float(args.d)
                and float(mydict[k][4]) >= float(args.l)):
            my_selected += 1
            Genus_slc=build_dict(k, Genus_slc)

            print("\t".join(mydict[k]), file=fout)

        print("\t".join(id_hit[k]), file=fout2)

    return Genus, Genus_slc, my_selected

def counter(genus_dict, fout):
    t_sp = 0 #total species
    genome_count = 0 #total genomes
    for g1, s1 in genus_dict.items():
        t_sp += len(s1.keys())
        for st1, ns in s1.items():
            n = len(ns.keys())
            print("## {} {} {}".format(g1, st1, n), file=fout)
            genome_count += n
    return t_sp, genome_count

def select_hits_and_printing_out_results(file_out, dict_id_hit):

    with open(file_out + "_selected", "w") as fout, open(file_out, "w") as fout2 :
        print("# {}".format("\t".join(header)), file=fout)
        print("# {}".format("\t".join(header)), file=fout2)
        Genus, Genus_slc, my_selected=get_hits_by_genus(dict_id_hit, fout, fout2)

        print("## Total number of unique best hits: {}".format(len(dict_id_hit.keys())), file=fout2)

        print("## Genus Sp strains", file=fout2)

        total_genus = len(Genus.keys())
        total_spp, genome_counter = counter(Genus, fout2)

        print("## Total number of genus: {} species :{} - strains: {}".format(total_genus, total_spp, genome_counter), file=fout2)
        print("       Total number of genus: {} species :{} - strains: {}".format(
                total_genus, total_spp, genome_counter))
        print("       ---------")
        print("## Total number of unique hits above threshold (% identity >= {} - % Query_cov >= {}): {}".format(args.d, args.l, my_selected), file=fout)
        print("## Genus Spp  hits", file=fout)

        t_spp, genome_counter_slc = counter(Genus_slc, fout)
        t_genus = len(Genus_slc.keys())

        print("## Total number of genus: {} - species :{} - strains: {}".format(t_genus, t_spp, genome_counter_slc), file=fout)

    return my_selected, t_genus,t_spp, genome_counter_slc

header = [ "query acc.", "subject acc.", "% identity", "query length",
"% query coverage", "alignment length", "mismatches", "gap opens",
"q. start", "q. end", "s. start", "s. end", "evalue", "bit score"]

id_hit=storing_blast_hits(args.inf)
print("       Total number of unique best hits: {}".format(len(id_hit.keys())))

id_hit_selected, to_genus,to_spp, genome_counter_slc = select_hits_and_printing_out_results(args.out, id_hit)

print(
    "       Total number of unique hits above threshold (% identity >= {} - % Query_cov >= {}): {}".format(
        args.d, args.l, id_hit_selected))
print("       Total number of genus: {} - species :{} - strains: {}".format(to_genus,
                                                                            to_spp, genome_counter_slc))

