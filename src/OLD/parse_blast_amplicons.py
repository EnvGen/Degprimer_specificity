#!/usr/bin/env python3

import os
import argparse
import re

usage = 'parse_blast_amplicons.py -i -p -r -d -l -M -m -o -a -k'
description = 'This program prints the predicted amplicons from primers'
#args.inf =  Results/BlastnEv_evalue15_strandboth_taskblastn_word_size11_max_target_seqs1000/FwL201_Rv570vsVibrioNtdb.blastn
#args.p = Primer/FwL201_Rv570.fasta

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument(
    '-i',
    dest='inf',
    help='input file .blastn',
    default="/root/Vibrio_BLAST/Results/BlastnEv_evalue15_strandboth_taskblastn_word_size11_max_target_seqs1000/FwL201_Rv576vsVibrioNtdb.blastn")
# required=True)
parser.add_argument('-p', dest='p', help='primers fasta file',
                    default="/root/Vibrio_BLAST/Primer/FwL201_Rv576.fasta")
# required=True)
parser.add_argument(
    '-c',
    dest='c',
    help='Input fasta (non-interleave) file, complete sequences ',
    default="/root/Vibrio_BLAST/Vibrio_Complete_genomes/all_Vibrio_Complete_genomes.fasta")
# required=True)
parser.add_argument(
    '-d', dest='d', help='identity threshold default 100', default=100)
parser.add_argument(
    '-l',
    dest='l',
    help='percentage query coverage threshold default 100',
    default=100)
parser.add_argument('-m', dest='m', help='min amplicon length', default=50)
parser.add_argument(
    '-M',
    dest='M',
    help='Max amplicon length default 1000',
    default=1000)
parser.add_argument('-o', dest='o', help='output file, hits',
                    default="TEST.txt")
parser.add_argument('-a', dest='a', help='output file .fna, amplicons',
                    default="TEST.fna")
parser.add_argument(
    '-k',
    dest='k',
    help='prefix output file .tsv, krona',
    default='TEST')

args = parser.parse_args()


def revcomp(line):
    old_chars = "ACGT"
    replace_chars = "TGCA"
    new = line.translate(line.maketrans(old_chars, replace_chars))[::-1]
    return new


primers_seq = {}
with open(args.p, "r") as fp:
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            clean = line[1:].replace(" ", "")
        else:
            primers_seq[clean] = line


id_hit = {}
with open(args.inf, "r") as fin:
    for line in fin:
        line = line.rstrip()
# Fields: query acc., subject acc., % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
        # primer_0_FwL201_forward
        # LS997867.1_Vibrio_cholerae_strain_NCTC_30_genome_assembly
        # 100.000 20      0       0      120      177560  177579  0.029   37.4
        if not line.startswith("#"):

            line = line.split()
           # id=line[1].split("_")[0] #it doesn't work since same strain may
           # have two different AccNumber because there are two Chromosomes
            loci = line[1]

        #    id = re.sub("^.*_Vibrio", "Vibrio", loci, count=1)
            id = loci.split("_")
            id = "_".join(id[2:])
            id = re.sub(".chromosome.*$", "", id, count=1)
            id = re.sub(".genome.*$", "", id, count=1)
            id = re.sub(".plasmid.*$", "", id, count=1)
            id = re.sub(".DNA.*$", "", id, count=1)
            id = re.sub(".complete.*$", "", id, count=1)
            id = re.sub("sp.", "sp", id, count=1)
            identity = float(line[2])
            align_len = float(line[3])
            pn = line[0]
            start = line[8]
            stop = line[9]

            primer_len = len(primers_seq[pn])

            p_align = round(align_len / primer_len * 100, 2)
            if p_align >= float(args.l) and identity >= float(args.d):
                if loci in id_hit:
                    if id in id_hit[loci]:
                        id_hit[loci][id].append((pn, start, stop))
                    else:
                        id_hit[loci] = {id: [(pn, start, stop)]}

                else:
                    id_hit[loci] = {id: [(pn, start, stop)]}

G_species = {}
strains = set()
al = []
total_amplicons = 0
selected_hits = {}
with open(args.o, "w") as fout:
    for l in id_hit:
        # print(l) NZ_CP042300.1_Vibrio_cholerae_O1_strain_AAS91_chromosome_2
        for i in id_hit[l]:
            #            print(i) Vibrio_cholerae_O1_strain_AAS91
            ps = id_hit[l][i]
    #        print(ps) [('primer_42_P74_326_forward', '1812666', '1812650'), ('primer_2_PR186_1166_reverse', '1811807', '1811826')]

            rev = set()
            forw = set()
            for u in ps:
                #    print(u) ('primer_42_P74_326_forward', '1812666', '1812650')
                if re.search("reverse", u[0]):
                    rev.add(u)
                else:
                    forw.add(u)
            for q in forw:
                for w in rev:
                    if q and w:
                        a = int(q[1])
                        b = int(q[2])
                        c = int(w[1])
#                        d=int(w[2])
                        if a < b:
                            direction = "+"
                            sta = a
                            sto = c
                        else:
                            direction = "-"
                            sta = c
                            sto = a

                        amplicon_length = sto - sta
                        if (amplicon_length >= int(args.m)
                                and amplicon_length <= int(args.M)):

                            # 3
                            # Vibrio_parahaemolyticus_strain_FDAARGOS_667
                            gns = i.split("_")[0]  # Vibrio
                            sps = i.split("_")[1]  # parahaemolyticus
                            # strain_FDAARGOS_667
                            strs = "_".join(i.split("_")[2:])

#                            print(gns, sps, strs)
                            if gns in G_species:
                                if sps in G_species[gns].keys():
                                    if strs in G_species[gns][sps].keys():
                                        G_species[gns][sps][strs] += 1
                                    else:
                                        G_species[gns][sps][strs] = 1
                                else:
                                    G_species[gns][sps] = {strs: 1}
                            else:
                                G_species[gns] = {sps: {strs: 1}}

# 3
#                            st = i.split("_")[1]
#                            st=st.split()[0]
#                            if st in species:
#                                species[st] += 1
#                            else:
#                                species[st] = 1
                            if gns == "Vibrio":
                                total_amplicons += 1
                                strains.add(i)
                                al.append(amplicon_length)
                            # print(l,direction, sta, sto)
                            # NC_005140.1_Vibrio_vulnificus_YJ016_chromosome_II
                            # - 1811807 1812666
                            if l in selected_hits:
                                selected_hits[l].append((direction, sta, sto))
                            else:
                                selected_hits[l] = [(direction, sta, sto)]

                            print("*** Strain: {} - Loci: {}".format(i, l), file=fout)
                            print("forward {} - reverse {} - amplicon length {} - direction {} - start {} - stop {}".format(q, w, amplicon_length, direction, sta, sto), file=fout)
    if len(al) > 0:
        print("#            From genus Vibrio: Total amplicons {} - species {} - strains {}".format(total_amplicons, len(G_species["Vibrio"].keys()), len(strains)), file=fout)
        print("#                Avg amplicon size {} Max {} Min {}".format(round(sum(al) / len(al), 1), max(al), min(al)), file=fout)
    else:
        print("#            From genus Vibrio: Total amplicons {} - species {} - strains {}".format(total_amplicons, 0, 0), file=fout)

if len(al) > 0:
    print(
        "#                From genus Vibrio: Total amplicons {} - species {} - strains {}".format(
            total_amplicons, len(
                G_species["Vibrio"].keys()), len(strains)))
    print("#                    Avg amplicon size {} Max {} Min {}".format(
        round(sum(al) / len(al), 1), max(al), min(al)))


with open(args.c, "r") as fc, open(args.a, "w") as fam, open(args.k + "_all_amplicons_krona.tsv", "w") as fk, open(args.k + "_unique_amplicons_krona.tsv", "w") as fku:
    if len(selected_hits.keys()) > 0:

        seq_uniq = {}
        copy = False
        for line in fc:
            line = line.rstrip()
            if line.startswith(">"):
                # print(line)
                # #>NZ_LS997868.1_Vibrio_cholerae_strain_NCTC_30_chromosome_2
                # complete_sequence
                h = line.split()[0][1:]
    # print(h) NZ_LS997868.1_Vibrio_cholerae_strain_NCTC_30_chromosome_2
                if h in selected_hits:
                    header = re.sub(".chromosome.*$", "", line, count=1)
                    header = re.sub(".genome.*$", "", header, count=1)
                    header = re.sub(".plasmid.*$", "", header, count=1)
                    header = re.sub(".DNA.*$", "", header, count=1)
                    header = re.sub(".complete.*$", "", header, count=1)
                    header = re.sub("sp.", "sp", header, count=1)
                    header += " [direction +]"
                    k = h
                    copy = True
            else:
                if copy:
                    for item in selected_hits[k]:
                        #                    ini=int(selected_hits[k][0][1])-1
                        #                    final=int(selected_hits[k][0][2])
                        ini = int(item[1]) - 1
                        final = int(item[2])
    #                    if selected_hits[k][0][0] == "+":
                        head = header + "[length=" + str(final - ini) + "]"
                        if item[0] == "+":
                            print(head, file=fam)
                            seqprint = line[ini:final]
                            print(seqprint, file=fam)
    # print(head) >NC_016628.1_Vibrio_furnissii_NCTC_11218 [direction
    # +][length=860]
                            strname = head.split("_")
                            strname = "_".join(strname[2:])
    #                            strname=re.sub("^.*_Vibrio", "Vibrio",head,count=1)
                            strname = strname.split()[0]
    #                            print(strname) Vibrio_furnissii_NCTC_11218
                            if seqprint in seq_uniq:
                                seq_uniq[seqprint].append((strname))
                            else:
                                seq_uniq[seqprint] = [(strname)]

                        else:
                            print(head, file=fam)
                            seqtoprint = revcomp(line[ini:final])
                            print(seqtoprint, file=fam)
                            strname = head.split("_")
                            strname = "_".join(strname[2:])
        #                        strname=re.sub("^.*_Vibrio", "Vibrio",head,count=1)
                            strname = strname.split()[0]
                            if seqtoprint in seq_uniq:
                                seq_uniq[seqtoprint].append((strname))
                            else:
                                seq_uniq[seqtoprint] = [(strname)]

                        copy = False

        print("#            Total number of unique sequences (100% identity, same size): {}".format(len(seq_uniq.keys())), file=fam)
        print("#                Total number of unique sequences (100% identity, same size): {}".format(
            len(seq_uniq.keys())))
    # for v in seq_uniq.values():
    #    print(v)

        identifiable_strain = set()
        for names in seq_uniq.values():
            #    print(names)
            if len(names) == 1:
                identifiable_strain.add(names[0])

        Genus = {}
        total_vibrio_identyf_sp = set()
        total_vibrio_identyf_strn = set()
        for sro in identifiable_strain:
                # Vibrio_coralliilyticus_strain_RE22
            gns = sro.split("_")[0]
            sps = sro.split("_")[1]
            strs = "_".join(sro.split("_")[2:])
            if gns == "Vibrio":
                total_vibrio_identyf_sp.add(sps)
                total_vibrio_identyf_strn.add(gns + "_" + sps + "_" + strs)
        #        print(gns, sps, strs)
        #        if s in species_hits:
        #            species_hits[s] += 1
        #        else:
        #            species_hits[s] = 1
            if gns in Genus:
                if sps in Genus[gns].keys():
                    if strs in Genus[gns][sps].keys():
                        Genus[gns][sps][strs] += 1
                    else:
                        Genus[gns][sps][strs] = 1
                else:
                    Genus[gns][sps] = {strs: 1}
            else:
                Genus[gns] = {sps: {strs: 1}}
        # print(Genus)
        if "Vibrio" in Genus.keys():
            print(
                "#                From genus Vibrio: Total number of strains 100% identifiable: {} strains from {} species".format(
                    len(total_vibrio_identyf_strn),
                    len(total_vibrio_identyf_sp)))
            print("#                    particularly,")

            if "cholerae" in Genus["Vibrio"].keys():
                print("#                                Cholerae {} strains".format(
                    len(Genus["Vibrio"]["cholerae"].keys())))
            else:
                print("#                                Cholerae 0 strain")

            if "vulnificus" in Genus["Vibrio"].keys():
                print("#                                vulnificus {} strains".format(
                    len(Genus["Vibrio"]["vulnificus"].keys())))
            else:
                print("#                                vulnificus 0 strain")

            if "parahaemolyticus" in Genus["Vibrio"].keys():
                print("#                                parahaemolyticus {} strains".format(
                    len(Genus["Vibrio"]["parahaemolyticus"].keys())))
            else:
                print("#                                parahaemolyticus 0 strain")

            if "alginolyticus" in Genus["Vibrio"].keys():
                print("#                                alginolyticus {} strains".format(
                    len(Genus["Vibrio"]["alginolyticus"].keys())))
            else:
                print("#                                alginolyticus 0 strain")

            if "sp" in Genus["Vibrio"].keys():
                print("#                                sp {} strains".format(
                    len(Genus["Vibrio"]["sp"].keys())))
            else:
                print("#                                sp 0 strain")
            print("Identifiable Vibrio strains: {}".format(
                ",".join(total_vibrio_identyf_strn)))

        print("         -------------------------------------")
        contador_g = 0
        contador_sp = 0
        contador_str = 0
        for generos in Genus:
            if generos != "Vibrio":
                contador_g += 1
                for speci in Genus[generos]:
                    contador_sp += 1
                    for strain in Genus[generos][speci]:
                        contador_str += 1
                        print("{}\t{}\t{}\t{} {} {}".format(Genus[generos][speci][strain], generos, speci, generos, speci, strain), file=fku)
            else:
                for speci in Genus[generos]:
                    for strain in Genus[generos][speci]:
                        print("{}\t{}\t{}\t{} {} {}".format(Genus[generos][speci][strain], generos, speci, generos, speci, strain), file=fku)

        print(
            "# 100% identifiable amplicons that are not Vibrio: Genus {} Spp {} strains {}".format(
                contador_g,
                contador_sp,
                contador_str))

        contador_g = 0
        contador_sp = 0
        contador_str = 0
        for generos in G_species:
            if generos != "Vibrio":
                contador_g += 1
                for speci in G_species[generos]:
                    contador_sp += 1
                    for strain in G_species[generos][speci]:
                        contador_str += 1
                        print("{}\t{}\t{}\t{} {} {}".format(G_species[generos][speci][strain], generos, speci, generos, speci, strain), file=fk)
            else:
                for speci in G_species[generos]:
                    for strain in G_species[generos][speci]:
                        print("{}\t{}\t{}\t{} {} {}".format(G_species[generos][speci][strain], generos, speci, generos, speci, strain), file=fk)

        print(
            "# Number of amplicons that are not Vibrio: Genus {} Spp {} strains {}".format(
                contador_g,
                contador_sp,
                contador_str))
        if contador_g > 0:
            print("    where, genus (exluding Vibrio) are: {}".format(
                ",".join([g for g in G_species.keys() if g != "Vibrio"])))

        print("         -------------------------------------")

    else:
        print("#            no amplicons finds or amplicon size longer/smaller than threshold ({},{})".format(args.m, args.M))
        print("#            no amplicons finds or amplicon size longer/smaller than threshold ({},{})".format(args.m, args.M), file=fam)
