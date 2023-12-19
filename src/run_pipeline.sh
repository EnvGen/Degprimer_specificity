#!/bin/bash -l
wkd=$(pwd)

err_report() {
    echo "Error on line $1 - script run_pipeline.sh"
    cd $wd
    mess2="Use the command:\n    'bash src/run_pipeline.sh unlock' "
    echo -e "TIP $mess2 and run the pipeline again"
    if test -f "tempo"; then rm tempo; fi
    exit 1
}
trap 'err_report $LINENO' ERR

configFile=$wkd/support_files/primers_config.json
#Converting json to bash configuration
if [ -s $configFile ]; then
cat $configFile | sed s/"[{|}]"//g | sed s/":"/"="/g | sed s/" "//g | sed s/",$"//g | sed s/'"'//g > tempo
. tempo

fi

####
#Checking key parameters setting
if [[ "$workdir" != /* ]]; then
    echo "Please provide an absoltute path to the working directory in configfile e.g., 'workdir': '/absolute/path/to/working_directory/' "
    echo "Please double-check the absolute path to working directory $workdir and to configfile $configFile"
    rm tempo
    exit 1
fi

options=("$output_dir_name" "$selected_genus" "$target" "$primer_list_file")
for o in "${options[@]}"; do if [ -z "$o" ]; then echo "A key parameter is undefined. Please check in the $configFile file the parameters used"; rm tempo; exit 1; fi; done

if  [ -z "$dir_database" ] && [ -z "$database" ] ; then
    echo "Please provide either the path to the directory containing the genomes.gz files or to the fasta file that includes all the genomes. Please check in the $configFile file the parameters used"; rm tempo; exit 1
fi

####


if [ "$1" != "unlock" ]; then
  bp=$( echo $blast_params | sed s/'\-'/'_'/g | sed s/"'"//g | sed s/'_+'/'_'/g )
  params_dir=$( echo "id_"$identity"_QC_"$Query_cov | sed s'/,//g' )$bp

  databasefasta=$database

  if [ ! -s $databasefasta ]; then
        echo "INFO: Extracting sequences from database"
        ls $dir_database/*genomic.fna.gz | while read f; do gunzip -cd $f >> $dir_database/tmpDB; done
        cat $dir_database/tmpDB | sed "s/ /_/g" | sed "s/,_/ /g" > $databasefasta
        rm $dir_database/tmpDB

  fi
# make sure the fatsa file is non interleaved
  check_interleaved=$(head -3 $databasefasta | grep -c ">")
  if [ "$check_interleaved" -ne 2 ]; then
        echo "INFO: Converting $databasefasta file to non-interleaved"
        awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $databasefasta > $wkd/"tmpDB_noninerleaved"
        mv $wkd/"tmpDB_noninerleaved" $databasefasta

  else
        echo "INFO: $databasefasta file is non-interleaved"

  fi

### Creating blast_db
  s=$(basename $database)
  blast_dbname="${s%.*}"

  if [ ! -d "$wkd/blast_db/$blast_dbname" ] || [ -z "$(ls -A $wkd/blast_db/$blast_dbname/)" ]; then
      echo "INFO: Creating blast dababase from $s"
       mkdir -p $wkd/blast_db/$blast_dbname
       makeblastdb -in "$databasefasta" -out "$wkd"/blast_db/"$blast_dbname"/Ntdb -dbtype nucl
  fi

  mkdir -p $wkd/$output_dir_name/$params_dir

  echo "# Primer specificity analysis" > $wkd/$output_dir_name/$params_dir/Summary.txt

fi


cat $primer_list_file | sed  "s/ \+/ /g" | while read line; do
if [ $(echo "$line" | grep "#" | wc -l) == 0 ] && [ $(echo "$line" | wc -l) != 0 ]; then
      id=$( echo "$line" | cut -d " " -f1 )
      namef=$( echo "$line" | cut -d " " -f2 )
      seqsf=$( echo "$line" | cut -d " " -f3 )
      namer=$( echo "$line" | cut -d " " -f4 )
      seqsr=$( echo "$line" | cut -d " " -f5 )

      if [ "$1" == "unlock" ]; then
        snakemake -s support_files/Primers.smk --config primer=$id seqf=$seqsf namef=$namef seqr=$seqsr namer=$namer --unlock
      else

        echo "INFO: Primer specificity analysis - Target gene: $id - forward primer: $namef $seqsf - reverse primer: $namer $seqsr"
        echo "#Target gene: $id - forward primer: $namef $seqsf - reverse primer: $namer $seqsr" >> $wkd/$output_dir_name/$params_dir/Summary.txt
        snakemake -s support_files/Primers.smk --cores $threads --config primer=$id seqf=$seqsf namef=$namef seqr=$seqsr namer=$namer --quiet >> $wkd/$output_dir_name/$params_dir/Summary.txt

      fi
fi
done


if [ "$1" != "unlock" ]; then
    echo "Printing out the best primers list: $wkd/$output_dir_name/$params_dir/SPECIAL_primers_list.txt"
    if  [ -z "$list_of_excl_genus" ]; then list_of_excl_genus="''"; fi
    python $wkd/src/parse_Summary.py -i $wkd/$output_dir_name/$params_dir/Summary.txt -o $wkd/$output_dir_name/$params_dir/SPECIAL_primers -d $target -s $min_idt_species -t $min_idt_strains -g $max_idt_genus -n $selected_genus -x $list_of_excl_genus
fi
rm tempo
