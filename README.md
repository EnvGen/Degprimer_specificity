The pipeline evaluates the specificity of degenerated primers -- Designed for genes from Vibrio genus, and support NCBI files

## Installation ##
Dowload the pipeline using the command:

        git clone https://github.com/lfdelzam/degprimer_specificity.git

The pipeline uses the programs:

[snakemake](https://snakemake.github.io) 3.13.3

[blastn](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs) v2.5.0

[krona](https://github.com/marbl/Krona) v2.7.1

## Create conda environment ##

conda create -n primers -c bioconda snakemake=3.13.3 blastn=2.5.0 krona=2.7.1

## Usage ##

set pipeline parameters in primers_config.json using the command:

    nano support_files/primers_config.json

and modify the parameters and save changes by taping `ctrl x` and tape `y`:

    "workdir": "/abs/path/to/degprimer_specificity",
    "dir_database": "/abs/path/to/Complete_genomes/", --Directory containing the NCBI genomes, no required if database file provided --
    "database": "/abs/path/to/complete_genomes.fasta", -- name required. if file not provided, it will be generated using dir_database --
    "output_dir_name": "your option",  -- required --
    "threads":12, -- used during blastn ---
    "blast_output_option":"-outfmt 7",
    "blast_params": "-evalue 20 -strand 'both' -task 'blastn' -word_size 11 -max_target_seqs 1000",
    "identity":"99",  -- The percent identity describes how similar the primer sequence is to the NCBI hit --
    "Query_cov": "99", -- the % of the primer length that aligns with the NCBI hit --
    "Max_ampl_size": "1100",
    "min_ampl_size": "200",
    "min_idt_species": 48, -- Used when selecting primers. Minimum number of species with unique amplicon that can be identified --
    "min_idt_strains": 153, -- Used when selecting primers. Minimum number of strains with unique amplicon that can be identified --
    "max_idt_genus": 50, -- Used when selecting primers. Maximum number of non target-genus from which amplicons are generated --
    "selected_genus": "your option" --targeted genus, e.g., "Vibrio" --
    "list_of_selected_spps": "your options" --List of specifique species you want to count, e.g., "cholerae,vulnificus,parahaemolyticus,alginolyticus,sp"--
    "target": "your option", -- level at which you want to check the primer specificity in the target database. If you want to check only on the selected genus type the genus name used in "selected_genus" (e.g., "Vibrio"). If you want to check on other genus, type "all" --
    "primer_list_file": "/abs/path/to/List_of_primers.txt" -- List of primers to be analysed --
