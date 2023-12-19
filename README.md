The pipeline evaluates the specificity of degenerated primers -- Designed for genes from a target genus, and support NCBI files

## Installation ##

Download the pipeline using the command:

        git clone https://github.com/lfdelzam/degprimer_specificity.git

The pipeline uses the programs:

[snakemake](https://snakemake.github.io) 3.13.3

[blastn](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs) v2.5.0

[krona](https://github.com/marbl/Krona) v2.7.1

## Create conda environment ##

conda create -n primers -c bioconda snakemake=3.13.3 blast=2.5.0 krona=2.7.1

## Usage ##

set pipeline parameters in primers_config.json using the command:

    nano support_files/primers_config.json

and modify the parameters and save changes by taping `ctrl x` and tape `y`:

    "workdir": "/abs/path/to/degprimer_specificity",
    "dir_database": "/abs/path/to/Complete_genomes/", --Directory containing the NCBI genomes, no required if database file provided --
    "database": "/abs/path/to/complete_genomes.fasta", -- name required. if file not provided, it will be generated using dir_database --
    "output_dir_name": "your option",  -- Name of the output directory, required --
    "threads":12, -- used during blastn ---
    "blast_output_option":"-outfmt 7",
    "blast_params": "-evalue 20 -strand 'both' -task 'blastn' -word_size 11 -max_target_seqs 1000",
    "identity":"99",  -- The percent identity describes how similar the primer sequence is to the NCBI hit --
    "Query_cov": "99", -- the % of the primer length that aligns with the NCBI hit --
    "Max_ampl_size": "1100",
    "min_ampl_size": "50",
    "min_idt_species": 1, -- Used when selecting primers. Minimum number of species with unique amplicon that can be identified --
    "min_idt_strains": 1, -- Used when selecting primers. Minimum number of strains with unique amplicon that can be identified --
    "max_idt_genus": 1, -- Used when selecting primers. Maximum number of non target-genus from which amplicons are generated --
    "selected_genus": "your option" --target-genus, e.g., "Vibrio" -- Required
    "list_of_special_spps": "your options" --List of specifique species you want to count, e.g., "cholerae,vulnificus,parahaemolyticus,alginolyticus,sp". Use "" if no list--
    "list_of_excl_genus": "your list", -- List of genus from which primers must not generate amplicons, e.g., "Aliivibrio". Use "" if no list --
    "target": "your option", -- level at which you want to check the primer specificity in the target database. If you want to check only on the selected genus type the genus name used in "selected_genus" (e.g., "Vibrio"). If you want to check on other genus, type "all". Required --
    "primer_list_file": "/abs/path/to/List_of_primers.txt" -- List of primers to be analysed. Required --


Activate environment and run the pipeline using the command:

    conda activate primers
    bash src/run_pipeline.sh


# Output

in <output_dir_name>:
├── BlastnEv_evalue<your_option>_strandboth_taskblastn_word_size<your_option>_max_target_seqs<your_option> -- blast table results for each primers pair evaluated
├── Spp_and_strains.txt -- list and number of Species/ strains present in the reference databse used to check the primer specificity
├── amplicons - Directory containing the amplicons sequences genereated for each primers pair evaluated
├── id_<your_option>_QC_<your_option>_evalue<your_option>_strandboth_taskblastn_word_size<your_option>_max_target_seqs<your_option>
     ├── SPECIAL_primers -- List of primers passing the selection filters
     ├── SPECIAL_primers_list.txt -- List of primers passing the selection filters in a type List_of_primers.txt format  
     ├── Summary.txt -- Log file containing key information. Example of output:

            # Primer specificity analysis
            #Target gene: wzc_1 - forward primer: Primer_99_at_42 GACGAAATCGATTTRGGC - reverse primer: Primer_184_at_142 AGGYGTTGAAAGHAGYG
                    total number of unique sequences from wzc_1_Primer_99_at_42 forward: 2  -- Primer degeneracy
                    total number of unique sequences from wzc_1_Primer_184_at_142 reverse: 12 -- Primer degeneracy
                   Total number of unique best hits: 415 -- Total blast hits before filtering by %identity and Query coverage
                   Total number of genus: 35 species :44 - strains: 415 -- Total blast hit on the selected genus before filtering by %identity and Query coverage
                   ---------
                   Total number of unique hits above threshold (% identity >= 99 - % Query_cov >= 99): 413 -- Total blast hits after filtering by %identity and Query coverage, before evaluating if an amplicon can be created
                   Total number of genus: 33 - species :42 - strains: 413 -- Total blast hit on the selected genus after filtering by %identity and Query coverage, before evaluating if an amplicon can be created
            #                From genus Vibrio: Total amplicons 365 - species 1 - strains 353 -- Total number of amplicons generated from the selected genus
            #                    Avg amplicon size 136.4 Max 736 Min 116
            #                Total number of unique sequences (100% identity, same size): 15290 -- among all the amplicons generated
            #                From genus Vibrio: Total number of strains 100% identifiable: 336 strains from 1 species -- Based on unique amplicons generated from the selected genus
            #                    particularly,
            #                                cholerae 0 strain
            #                                vulnificus 336 strains
            #                                parahaemolyticus 0 strain
            #                                alginolyticus 0 strain
            #                                sp 0 strain
            Identifiable Vibrio strains: -- List of strains, based on unique amplicons generated from the selected genus
            -------------------------------------
            # 100% identifiable amplicons that are not Vibrio: Genus 0 Spp 0 strains 0 -- Based on unique amplicons generated from the other genra
            # Number of amplicons that are not Vibrio: Genus 0 Spp 0 strains 0 -- Total number of amplicons generated from the other genra
            -------------------------------------

     ├── <primer_number>_<primers_pair>_hits -- Intermediare file, -- blast table results for each primers pair evaluated
     ├── <primer_number>_<primers_pair>_hits_selected -- Blast table hits above  %identity and Query coverage threshold for each primers pair evaluated
└── krona -- Directory containing Krona html files, of amplicon/hits for each primers pair evaluated
