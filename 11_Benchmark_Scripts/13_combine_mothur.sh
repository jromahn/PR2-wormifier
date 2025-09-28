input_path="PhytoArk_Euka02_2024_MUC_results/"
output_path="12_MOTHUR_ASSIGNMENTS"

cd $input_path
 
cat 8_PhytoArk_Euka02_2024_final_community_MUC__mothur_sequences__1.11_Euka02_Taxonomy_DB1.wang.taxonomy \
     8_PhytoArk_Euka02_2024_final_community_MUC__mothur_sequences__2.11_Euka02_Taxonomy_DB1.wang.taxonomy \
     > ../$output_path/9_PhytoArk_Euka02_MUC__Taxonomy_amplicons__DB1.wang.taxonomy

cat 8_PhytoArk_Euka02_2024_final_community_MUC__mothur_sequences__1.11_Euka02_Taxonomy_DB2.wang.taxonomy \
     8_PhytoArk_Euka02_2024_final_community_MUC__mothur_sequences__2.11_Euka02_Taxonomy_DB2.wang.taxonomy  \
     > ../$output_path/9_PhytoArk_Euka02_MUC__Taxonomy_amplicons__DB2.wang.taxonomy

cat 8_PhytoArk_Euka02_2024_final_community_MUC__mothur_sequences__1.11_Euka02_Taxonomy_DB3.wang.taxonomy \
     8_PhytoArk_Euka02_2024_final_community_MUC__mothur_sequences__2.11_Euka02_Taxonomy_DB3.wang.taxonomy \
     > ../$output_path/9_PhytoArk_Euka02_MUC__Taxonomy_amplicons__DB3.wang.taxonomy

cat 8_PhytoArk_Euka02_2024_final_community_MUC__mothur_sequences__1.11_Euka02_Taxonomy_DB4.wang.taxonomy \
     8_PhytoArk_Euka02_2024_final_community_MUC__mothur_sequences__2.11_Euka02_Taxonomy_DB4.wang.taxonomy \
     > ../$output_path/9_PhytoArk_Euka02_MUC__Taxonomy_amplicons__DB4.wang.taxonomy

cd ..