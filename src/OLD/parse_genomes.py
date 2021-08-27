#!/usr/bin/env python3

import os
import argparse
import re

usage = 'parse_vibrio_genomes -i -o'
description = 'Retrive the number of unique strains and the number of strain for each specie from complete genomes database'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument('-i', dest='inf', help='input file', required=True)
parser.add_argument('-o', dest='out', help='output file', required=True)
args = parser.parse_args()

id_hit = {}
with open(args.inf, "r") as fin:
    for line in fin:
        line = line.rstrip()
        if line.startswith(">"):
            id = re.sub(".chromosome.*$", "", line, count=1)
            id = re.sub(".genome.*$", "", id, count=1)
            id = re.sub(".plasmid.*$", "", id, count=1)
            id = re.sub(".DNA.*$", "", id, count=1)
            id = re.sub(".complete.*$", "", id, count=1)
            id = re.sub("sp.", "sp", id, count=1)
            # print(id)
            if id in id_hit:
                id_hit[id] += 1
            else:
                id_hit[id] = 1

with open(args.out, "w") as fout:
    print("# Total number of sequences: {}".format(len(id_hit.keys())), file=fout)

#--- estoy trabajando aca
    G_species = {}
    for s, n in id_hit.items():
        #        print("# {} {}".format(s,n), file=fout)
        gns = s.split("_")[2]
        sp = s.split("_")[3]
        strn = "_".join(s.split("_")[4:])
#        print(s, gns, sp, strn)
####
        if gns in G_species:
            if sp in G_species[gns].keys():
                if strn in G_species[gns][sp].keys():
                    G_species[gns][sp][strn] += 1
                else:
                    G_species[gns][sp][strn] = 1
            else:
                G_species[gns][sp] = {strn: 1}
        else:
            G_species[gns] = {sp: {strn: 1}}

####
#    print("*****************", file=fout)
    print("# Vibrio_spp  strains", file=fout)

#    species={k: v for k, v in sorted(species.items(), key=lambda item: item[1], reverse=True)}

###
    contador_g = 0
    contador_sp = 0
    contador_str = 0
    contador_Vsp = 0
    contador_Vstr = 0
    for generos in G_species:
        if generos != "Vibrio":
            contador_g += 1
            for speci in G_species[generos]:
                contador_sp += 1
                for strain in G_species[generos][speci]:
                    contador_str += len(G_species[generos][speci].keys())
        else:
            for speci in G_species["Vibrio"]:
                contador_Vsp += 1
                print("  {} {}".format(speci, len(G_species["Vibrio"][speci].keys())), file=fout)
                contador_Vstr += len(G_species["Vibrio"][speci].keys())
###
    print("# Total number of Vibrio species :{} - strains: {}".format(contador_Vsp, contador_Vstr), file=fout)
    print("# Total number of sequences that are not Vibrio: Genus {} Spp {} strains {}".format(contador_g, contador_sp, contador_str), file=fout)

# print(G_species)
