path="11_Benchmark_Database_Versions"


primer_forward="TTTGTCTGSTTAATTSCG"
primer_reverse="CACAGACCTGTTATTGC"

length_filter=500
percentage_id=0.85

conda_env="CRABS_database"


date_log=$(date '+%Y-%m-%d')

readme_file="README_CRABS_run_${date_log}.txt"

source /opt/anaconda3/etc/profile.d/conda.sh               #make sure to be able to switch between conda envs
conda activate $conda_env

echo "Day: $date_log" > $readme_file
echo "PRIMER Forward: $primer_forward" >> $readme_file
echo "PRIMER Reverse: $primer_reverse" >> $readme_file

echo "Setting for the pairwaise global alignment " >> $readme_file
echo "Length: $length_filter" >> $readme_file
echo "Percentage ID: $percentage_id" >> $readme_file
echo "Speed: medium" >> $readme_file
echo "Filter method: strict" >> $readme_file
echo "Coverage: 0.95" >> $readme_file


## run crabs
for file in $( ls  $path/*_fake_header.fasta); do
    file_name=$( basename $file .fasta)
    file_name=$(echo "$file_name" | sed "s/Sequence/Amplicon/g")
    echo $file_name
    echo $file

    echo "" >> $readme_file
    echo "Input: $file" >> $readme_file
    

    ## insilico PCR
    crabs insilico_pcr --input $file --output "$path/${file_name}.fasta" --fwd $primer_forward --rev $primer_reverse --error 4.5 --t 10
    amount_insilico1=$( grep ">" "$path/${file_name}.fasta"  | wc -l )

    # filter out too  long sequences
    echo "Remove amplicons which are too long: $prefix "
    vsearch -sortbylength "$path/${file_name}.fasta"  --maxseqlength $length_filter --output "$path/${file_name}_L$length_filter.fasta"
    amount_insilico2=$( grep ">" "$path/${file_name}_L$length_filter.fasta" | wc -l )

    #regain
    echo "regain fitting amplicons: $prefix "
    crabs pga --input $file --output "$path/${file_name}_L${length_filter}_P${percentage_id}.fasta" --database "$path/${file_name}_L$length_filter.fasta" \
            --fwd $primer_forward --rev $primer_reverse --speed medium --percid $percentage_id --coverage 0.95 --filter_method strict
    amount_insilico3=$( grep ">" "$path/${file_name}_L${length_filter}_P${percentage_id}.fasta" | wc -l )

    printf "insilico PCR  : "$path/${file_name}.fasta" \t $amount_insilico1 sequences \n" >> $readme_file
    printf "insilico PCR  - length filtered: "$path/${file_name}_L$length_filter.fasta" \t$amount_insilico2 sequences \n" >> $readme_file
    printf "insilico PCR - regained: "$path/${file_name}_L${length_filter}_P${percentage_id}.fasta" \t $amount_insilico3 sequences \n" >> $readme_file
done

conda deactivate