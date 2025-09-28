#########################################
# Manual to assign obtools3 data in mothur
## 1.)  transforms obitools table output to necceassary counttable: execute 02_prepare_data_for_mothur.R ( with neccessary changed) 
#### "  and : has to be removed
## 2.) prepare count table via: sed 's|"||g' 01_obitools_results/Dino_JR_1_counttable.tsv | sed 's|:|_|g' >> 02_mothur_results/Dino_JR_1_counttable_mothur.tsv 
## 3.) copy files because 02_mothus_results will be workdir: cp 01_obitools_results/*.fasta 02_mothur_results/
## 4.) execute this mothur script via: mothur "13_README_mothur.txt"
######### MOTHUR Script
## note: no - in names allowed
set.dir(input=/PATH/TO/YOURDATA/)
# classify Euka02 db1 
#classify.seqs(fasta=PhytoArk_Euka02_2024_MUC_results/8_PhytoArk_Euka02_2024_final_community_MUC__mothur_sequences__1.fasta, count=PhytoArk_Euka02_2024_MUC_results/8_PhytoArk_Euka02_2024_final_community_MUC__mothur_counttable__1.tsv , reference=11_Benchmark_Database_Versions/11_Euka02_database_V1.fasta ,taxonomy=11_Benchmark_Database_Versions/11_Euka02_Taxonomy_DB1.tax , cutoff =80)
#classify.seqs(fasta=PhytoArk_Euka02_2024_MUC_results/8_PhytoArk_Euka02_2024_final_community_MUC__mothur_sequences__2.fasta, count=PhytoArk_Euka02_2024_MUC_results/8_PhytoArk_Euka02_2024_final_community_MUC__mothur_counttable__2.tsv , reference=11_Benchmark_Database_Versions/11_Euka02_database_V1.fasta ,taxonomy=11_Benchmark_Database_Versions/11_Euka02_Taxonomy_DB1.tax , cutoff =80)

# classify Euka02 db2
#classify.seqs(fasta=PhytoArk_Euka02_2024_MUC_results/8_PhytoArk_Euka02_2024_final_community_MUC__mothur_sequences__1.fasta, count=PhytoArk_Euka02_2024_MUC_results/8_PhytoArk_Euka02_2024_final_community_MUC__mothur_counttable__1.tsv , reference=11_Benchmark_Database_Versions/11_Euka02_database_V2.fasta ,taxonomy=11_Benchmark_Database_Versions/11_Euka02_Taxonomy_DB2.tax , cutoff =80)
#classify.seqs(fasta=PhytoArk_Euka02_2024_MUC_results/8_PhytoArk_Euka02_2024_final_community_MUC__mothur_sequences__2.fasta, count=PhytoArk_Euka02_2024_MUC_results/8_PhytoArk_Euka02_2024_final_community_MUC__mothur_counttable__2.tsv , reference=11_Benchmark_Database_Versions/11_Euka02_database_V2.fasta ,taxonomy=11_Benchmark_Database_Versions/11_Euka02_Taxonomy_DB2.tax , cutoff =80)
# classify Euka02 db3
#classify.seqs(fasta=PhytoArk_Euka02_2024_MUC_results/8_PhytoArk_Euka02_2024_final_community_MUC__mothur_sequences__1.fasta, count=PhytoArk_Euka02_2024_MUC_results/8_PhytoArk_Euka02_2024_final_community_MUC__mothur_counttable__1.tsv , reference=11_Benchmark_Database_Versions/11_Euka02_database_V3.fasta ,taxonomy=11_Benchmark_Database_Versions/11_Euka02_Taxonomy_DB3.tax , cutoff =80)
#classify.seqs(fasta=PhytoArk_Euka02_2024_MUC_results/8_PhytoArk_Euka02_2024_final_community_MUC__mothur_sequences__2.fasta, count=PhytoArk_Euka02_2024_MUC_results/8_PhytoArk_Euka02_2024_final_community_MUC__mothur_counttable__2.tsv , reference=11_Benchmark_Database_Versions/11_Euka02_database_V3.fasta ,taxonomy=11_Benchmark_Database_Versions/11_Euka02_Taxonomy_DB3.tax , cutoff =80)
# classify Euka02 db2
classify.seqs(fasta=PhytoArk_Euka02_2024_MUC_results/8_PhytoArk_Euka02_2024_final_community_MUC__mothur_sequences__1.fasta, count=PhytoArk_Euka02_2024_MUC_results/8_PhytoArk_Euka02_2024_final_community_MUC__mothur_counttable__1.tsv , reference=11_Benchmark_Database_Versions/11_Euka02_database_V4.fasta ,taxonomy=11_Benchmark_Database_Versions/11_Euka02_Taxonomy_DB4.tax , cutoff =80)
classify.seqs(fasta=PhytoArk_Euka02_2024_MUC_results/8_PhytoArk_Euka02_2024_final_community_MUC__mothur_sequences__2.fasta, count=PhytoArk_Euka02_2024_MUC_results/8_PhytoArk_Euka02_2024_final_community_MUC__mothur_counttable__2.tsv , reference=11_Benchmark_Database_Versions/11_Euka02_database_V4.fasta ,taxonomy=11_Benchmark_Database_Versions/11_Euka02_Taxonomy_DB4.tax , cutoff =80)
