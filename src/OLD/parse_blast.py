#!/usr/bin/env python3

import os
import argparse
import re

usage = 'parse_blast.py -i -n -f -r -o'
description = 'This program count the number of unique hits from a blast table, above an user-defined threshold'

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

header = [ "query acc.", "subject acc.", "% identity", "query length",
"% query coverage", "alignment length", "mismatches", "gap opens",
"q. start", "q. end", "s. start", "s. end", "evalue", "bit score"]

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

def store_max_min(k, id_hit):
    if float(id_hit[k][13]) < lb:
        lb = float(id_hit[k][13])
    if float(id_hit[k][12]) < le:
        le = float(id_hit[k][12])
    if float(id_hit[k][2]) < li:
        li = float(id_hit[k][2])
    if float(id_hit[k][4]) < lpa:
        lpa = float(id_hit[k][4])
    if float(id_hit[k][13]) > Mb:
        Mb = float(id_hit[k][13])
    if float(id_hit[k][12]) > Me:
        Me = float(id_hit[k][12])
    if float(id_hit[k][2]) > Mi:
        Mi = float(id_hit[k][2])
    if float(id_hit[k][4]) > Mpa:
        Mpa = float(id_hit[k][4])

    return lb, le, li, lpa, Mb, Me, Mi, Mpa                



id_hit=storing_blast_hits(args.inf)
print("       Total number of unique best hits: {}".format(len(id_hit.keys())))

id_hit_selected = {}
lowestevalue = le = 100
lowestbitscore = lb = 1000
lowestidentity = li = 101
lowestpa = lpa = 101
Maxevalue = Me = 0
Maxbitscore = Mb = 0
Maxidentity = Mi = 0
Maxpa = Mpa = 0
Genus = {}
Genus_slc = {}

with open(args.out + "_selected", "w") as fout, open(args.out, "w") as fout2, open("Summary_table", "a") as filetable:
    print("# {}".format("\t".join(header)), file=fout)
    print("# {}".format("\t".join(header)), file=fout2)

    for k in id_hit.keys():
        gns = k.split("_")[0]
        sps = k.split("_")[1]
        str = "_".join(k.split("_")[2:])
#        sps=sps.split()[0]
        if gns in Genus:
            if sps in Genus[gns].keys():
                if str in Genus[gns][sps].keys():
                    Genus[gns][sps][str] += 1
                else:
                    Genus[gns][sps][str] = 1
            else:
                Genus[gns][sps] = {str: 1}
        else:
            Genus[gns] = {sps: {str: 1}}

        if float(id_hit[k][13]) < lb:
            lb = float(id_hit[k][13])
        if float(id_hit[k][12]) < le:
            le = float(id_hit[k][12])
        if float(id_hit[k][2]) < li:
            li = float(id_hit[k][2])
        if float(id_hit[k][4]) < lpa:
            lpa = float(id_hit[k][4])

        if float(id_hit[k][13]) > Mb:
            Mb = float(id_hit[k][13])
        if float(id_hit[k][12]) > Me:
            Me = float(id_hit[k][12])
        if float(id_hit[k][2]) > Mi:
            Mi = float(id_hit[k][2])
        if float(id_hit[k][4]) > Mpa:
            Mpa = float(id_hit[k][4])
#        print(id_hit[k])
        if (float(id_hit[k][2]) >= float(args.d)
                and float(id_hit[k][4]) >= float(args.l)):
            id_hit_selected[k] = id_hit[k]
#            print(id_hit[k])
            gns_s = k.split("_")[0]
            sps_s = k.split("_")[1]
    #        sps=sps.split()[0]
            str_s = "_".join(k.split("_")[2:])
    #        sps=sps.split()[0]
            if gns_s in Genus_slc:
                if sps_s in Genus_slc[gns_s].keys():
                    if str_s in Genus_slc[gns_s][sps_s].keys():
                        Genus_slc[gns_s][sps_s][str_s] += 1
                    else:
                        Genus_slc[gns_s][sps_s][str_s] = 1
                else:
                    Genus_slc[gns_s][sps_s] = {str_s: 1}
            else:
                Genus_slc[gns_s] = {sps_s: {str_s: 1}}

            print("\t".join(id_hit[k]), file=fout)
            if float(id_hit[k][13]) < lowestbitscore:
                lowestbitscore = float(id_hit[k][13])
            if float(id_hit[k][12]) < lowestevalue:
                lowestevalue = float(id_hit[k][12])
            if float(id_hit[k][2]) < lowestidentity:
                lowestidentity = float(id_hit[k][2])
            if float(id_hit[k][4]) < lowestpa:
                lowestpa = float(id_hit[k][4])

            if float(id_hit[k][13]) > Maxbitscore:
                Maxbitscore = float(id_hit[k][13])
            if float(id_hit[k][12]) > Maxevalue:
                Maxevalue = float(id_hit[k][12])
            if float(id_hit[k][2]) > Maxidentity:
                Maxidentity = float(id_hit[k][2])
            if float(id_hit[k][4]) > Maxpa:
                Maxpa = float(id_hit[k][4])

        print("\t".join(id_hit[k]), file=fout2)

    print("## Total number of unique best hits: {}".format(len(id_hit.keys())), file=fout2)
    if len(Genus.values()) == 0:
        le = Me = " "
        lb = li = lpa = Mpa = 0

    print("##   Bitscore range ({},{})- Evalue range ({}, {}) - % identity range ({},{}) - % Query_cov range ({},{})".format(lb, Mb, le, Me, li, Mi, lpa, Mpa), file=fout2)
    print(
        "       Bitscore range ({},{})- Evalue range ({}, {}) - % identity range ({},{}) - % Query_cov range ({},{})".format(
            lb,
            Mb,
            le,
            Me,
            li,
            Mi,
            lpa,
            Mpa))

#    tempo_names=args.inf.split("vs")[0]
#    tempo_names=tempo_names.split("/")[-1]
#    namef=tempo_names.split("_")[0]
#    namer=tempo_names.split("_")[1]

#    Genus={k: v for k, v in sorted(Genus.items(), key=lambda item: item[1].values(), reverse=True)}
#    Genus_slc={k: v for k, v in sorted(Genus_slc.items(), key=lambda item: item[1].values(), reverse=True)}

    print("## Genus Sp strains", file=fout2)
    genome_counter = 0

    total_spp = 0

    total_genus = len(Genus.keys())
    for g, s in Genus.items():
        total_spp += len(s.keys())

        for st, ns in s.items():
            n = len(ns.keys())
#            print(g,st,n)
            print("## {} {} {}".format(g, st, n), file=fout2)
            genome_counter += n

    print("## Total number of genus: {} species :{} - strains: {}".format(total_genus, total_spp, genome_counter), file=fout2)
    print(
        "       Total number of genus: {} species :{} - strains: {}".format(
            total_genus,
            total_spp,
            genome_counter))
    print("       ---------")

    print("## Total number of unique hits above threshold (% identity >= {} - % Query_cov >= {}): {}".format(args.d, args.l, len(id_hit_selected.keys())), file=fout)
    print("##   Bitscore range ({},{})- Evalue range ({}, {}) - % identity range ({},{}) - % Query_cov range ({},{})".format(lowestbitscore, Maxbitscore, lowestevalue, Maxevalue, lowestidentity, Maxidentity, lowestpa, Maxpa), file=fout)
    print("## Genus Spp  hits", file=fout)
    genome_counter_slc = 0
    if len(Genus_slc.values()) == 0:
        lowestevalue = Maxevalue = " "
        lowestbitscore = lowestidentity = lowestpa = Maxpa = 0

    t_spp = 0
    t_genus = len(Genus_slc.keys())
    for g, s in Genus_slc.items():
        t_spp += len(s.keys())
        for st, ns in s.items():
            n = len(ns.keys())
            print("## {} {} {}".format(g, st, n), file=fout)
            genome_counter_slc += n

    print("## Total number of genus: {} - species :{} - strains: {}".format(t_genus, t_spp, genome_counter_slc), file=fout)

# print("{}\t{}\t{}\tTotal number of genus: {} - species :{} - strains:
# {}\tTotal number of genus: {} - species :{} - strains:
# {}".format(args.n, namef, namer, total_genus,total_spp, genome_counter,
# t_genus,t_spp, genome_counter_slc), file=filetable)
    print("{}\t{}\t{}\tTotal number of genus: {} - species :{} - strains: {}\tTotal number of genus: {} - species :{} - strains: {}".format(args.n, args.nf, args.nr, total_genus, total_spp, genome_counter, t_genus, t_spp, genome_counter_slc), file=filetable)
    print("{}\t{}\t{}\tBitscore range ({},{})- Evalue range ({}, {}) - % identity range ({},{}) - % Query_cov range ({},{})\tBitscore range ({},{})- Evalue range ({}, {}) - % identity range ({},{}) - % Query_cov range ({},{})".format("", "", "", lb, Mb, le, Me, li, Mi, lpa, Mpa, lowestbitscore, Maxbitscore, lowestevalue, Maxevalue, lowestidentity, Maxidentity, lowestpa, Maxpa), file=filetable)

print(
    "       Total number of unique hits above threshold (% identity >= {} - % Query_cov >= {}): {}".format(
        args.d, args.l, len(
            id_hit_selected.keys())))
print("       Bitscore range ({},{})- Evalue range ({}, {}) - % identity range ({},{}) - % Query_cov range ({},{})".format(
    lowestbitscore, Maxbitscore, lowestevalue, Maxevalue, lowestidentity, Maxidentity, lowestpa, Maxpa))
print("       Total number of genus: {} - species :{} - strains: {}".format(t_genus,
                                                                            t_spp, genome_counter_slc))
