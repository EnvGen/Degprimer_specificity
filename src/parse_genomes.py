#!/usr/bin/env python3

import os
import argparse
import re

usage = 'parse_genomes.py -i -o -g'
description = 'It retrives the number of unique strains and the number of strain for each specie from complete genomes database NCBI format'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument('-i', dest='inf', help='input file', required=True)
parser.add_argument('-o', dest='out', help='output file', required=True)
parser.add_argument('-g', dest='g', help='selected genus', default="Vibrio")

args = parser.parse_args()

def get_hits(file):
    hits = set()
    with open(file, "r") as fin:
        for line in fin:
            line = line.rstrip()
            if line.startswith(">"):
                id = re.sub(".chromosome.*$", "", line, count=1)
                id = re.sub(".genome.*$", "", id, count=1)
                id = re.sub(".plasmid.*$", "", id, count=1)
                id = re.sub(".DNA.*$", "", id, count=1)
                id = re.sub(".complete.*$", "", id, count=1)
                id = re.sub("sp.", "sp", id, count=1)
                hits.add(id)
    return hits

def Genome_counter(hits):
    Genus_species_strains = {}
    for s in hits:
        gns = s.split("_")[2]
        sp = s.split("_")[3]
        strn = "_".join(s.split("_")[4:])

        if gns in Genus_species_strains:
            if sp in Genus_species_strains[gns].keys():
                if strn in Genus_species_strains[gns][sp].keys():
                    Genus_species_strains[gns][sp][strn] += 1
                else:
                    Genus_species_strains[gns][sp][strn] = 1
            else:
                Genus_species_strains[gns][sp] = {strn: 1}
        else:
            Genus_species_strains[gns] = {sp: {strn: 1}}

    return Genus_species_strains

def print_out(file_out, hits, gen_counts, selected_genus):

    with open(file_out, "w") as fout:
        print("# Total number of sequences: {}".format(len(hits)), file=fout)
        print("# {}_spp  strains".format(selected_genus), file=fout)

        contador_g = 0
        contador_sp = 0
        contador_str = 0
        contador_Vsp = 0
        contador_Vstr = 0
        for generos in gen_counts:
            if generos != selected_genus:
                contador_g += 1
                for speci in gen_counts[generos]:
                    contador_sp += 1
                    for strain in gen_counts[generos][speci]:
                        contador_str += len(gen_counts[generos][speci].keys())
            else:
                for speci in gen_counts[selected_genus]:
                    contador_Vsp += 1
                    print("  {} {}".format(speci, len(gen_counts[selected_genus][speci].keys())), file=fout)
                    contador_Vstr += len(gen_counts[selected_genus][speci].keys())

        print("# Total number of {} species :{} - strains: {}".format(selected_genus, contador_Vsp, contador_Vstr), file=fout)
        print("# Total number of sequences that are not {}: Genus {} Spp {} strains {}".format(selected_genus, contador_g, contador_sp, contador_str), file=fout)


id_hit=get_hits(args.inf)
G_species=Genome_counter(id_hit)
print_out(args.out, id_hit, G_species, args.g)
