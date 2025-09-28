

#######################
# Aim: download species from NBCI to create a Reference Database, if species is not represented find sequences from same genus
# have to be added check if fasta file is empty
#######################

import os
import sys
import re
import numpy
from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
from datetime import datetime
import shutil
import time
from io import StringIO


# Set a true global timeout for all Entrez calls
import urllib
# Correct way: return raw binary stream — Entrez handles parsing
def _patched_open(*args, **kwargs):
    kwargs['timeout'] = 30  # or your TIMEOUT
    resp = urllib.request.urlopen(*args, **kwargs)
    content_type = resp.info().get_content_charset() or "utf-8"
    # wrap as text if retmode=text
    return io.TextIOWrapper(resp, encoding=content_type)




# define NCBI timeout time 
TIMEOUT=30
CONFIG_file="00_login_data.ini"


## read in config file with email
import configparser
config = configparser.ConfigParser()
with open(CONFIG_file) as f:
    config.read_string("[DEFAULT]\n" + f.read())

Entrez.email = config["DEFAULT"]["Entrez_email"]

###### extract date
date = datetime.now()
date = date.strftime("%d_%b_%Y")
date = str(date)


### Steps
FLAG_search_spec=True ### Step 1 - if species list should be 
FLAG_search_gen=True ### Step 2 - check if other species of the same genus are represented (depends on Step 1)

#ncbi search parameter 
search=" AND (18S OR small ribosomal) NOT (chloroplast[filter] OR mitochondrion[filter]) AND (350[SLEN] : 5000[SLEN])"
ncbi_db="nucleotide"

#input files are:genus,species_name\n (have to contain a header)
Input_folder="01_intermediate_results/"
file_missing = Input_folder+"4.1_Missing_Species.csv" ## species missing in database (Step 1)
file_missing_but_genus=Input_folder+"4.2_Present_Species.csv"  ## complete species list in pr2 for genera which still miss species(step 1 & 2)
file_missing_but_genus_new = Input_folder+"4.2_Present_Species_completed_withNCBI.csv" 
shutil.copy(file_missing_but_genus, file_missing_but_genus_new )


#output
NCBI_folder="02_NCBI_"+ date +"/" # where to store the downloaded files
NCBI_results_folder="02_NCBI_"+ date +"_results/" # where to store the downloaded files
output_prefix="05_Missing_PR2" # prefix of the output files

#######################################################################################################################
################################### PLEASE DON'T CHANGE ANYTHING FROM HERE ON ########################################
#######################################################################################################################

################ VARIABLES
#output
file_still_missing= NCBI_results_folder + output_prefix +"_Species_and_missing_NCBI__"+ date +".csv"
file_missing_genera= NCBI_results_folder +output_prefix +"_incomplete_genera_despite_NCBI__"+ date +".csv"
file_missing_genera2= NCBI_results_folder +output_prefix +"_MISSING_genera_despite_NCBI__"+ date +".csv"
file_asseccion_info= NCBI_results_folder +output_prefix +"_Downloaded_species_NCBI_info__"+ date +".tsv"
file_LOG= NCBI_results_folder + output_prefix + "LOG__" + date +".txt"
file_LOG2=NCBI_results_folder + output_prefix + "LOG__" + date +"_successfull_genera.txt"

if not os.path.exists(NCBI_folder):
   os.makedirs(NCBI_folder)

if not os.path.exists(NCBI_results_folder):
   os.makedirs(NCBI_results_folder)
################################################################################################################

# set so Entrez is underlying the try if there is long response or similar
#Entrez._open = _patched_open

#1.) check if species exists in ncbi and download fasta
if FLAG_search_spec:
    with open(file_missing, "r") as f, open(file_still_missing, "w") as f_new, open(file_missing_but_genus_new, "a") as f_new2,open(file_LOG, "a") as f_log, open(file_asseccion_info, "w") as f_assec:
        lines = f.readlines()[1:]
        for line in lines:
                line = line.strip() # remove whtie space and newline character from beginning and the end
                line = line.replace("\"","")
                line_information = line.split(",")
                print(line_information[1])
                new_search = line_information[1]+"[ORGANISM]"+ search
                f_log.write(new_search+ "\n") 
                attempts = 0
                while attempts < 3:
                    try:
                        print(attempts+1)
                        # check if species exist on ncbi
                        handle = Entrez.esearch(db=ncbi_db, retmax=10,term= new_search, timeout=TIMEOUT)
                        record = Entrez.read(handle)
                        handle.close()

                        # check if ncbi request was positive, if not continue
                        if len(record['IdList']) == 0:
                            f_new.write(line+ "\n") 
                            break

                        print("try to download... " )
                        id_list = list(record['IdList'])
                        out_file = NCBI_folder + line_information[1].replace(" ", "_") + "__NCBI.fasta"

                        # Download FASTA, should now catch the try function 
                        with Entrez.efetch(db=ncbi_db, id=id_list, rettype="fasta", retmode="text") as handle_fasta:
                            fasta_text = handle_fasta.read()  # Already str, no need to decode
                            all_records = SeqIO.parse(StringIO(fasta_text), 'fasta')
                            count = SeqIO.write(all_records, out_file, "fasta-2line")

                        # Check if file was created
                        if os.path.isfile(out_file):
                            f_new2.write(line_information[0] + "," + line_information[1] + "\n")

                            # Download GenBank and extract taxonomy info
                            with Entrez.efetch(db=ncbi_db, id=id_list, rettype="genbank", retmode="text") as handle_genbank:
                                genbank_text = handle_genbank.read()  # Already str, no need to decode
                                for record in SeqIO.parse(StringIO(genbank_text), 'genbank'):
                                    f_assec.write(
                                        line_information[1] + "\t" +
                                        record.annotations.get('organism', 'N/A') + "\t" +
                                        record.id + "\t" +
                                        out_file + "\t" +
                                        ";".join(record.annotations.get('taxonomy', [])) + "\n"
                                    )

                        break # success
                    except Exception as e:
                        print(f"Something went wrong with line '{line.strip()}': {e}")
                        f_log.write(line + f"\tproblem: {type(e).__name__} - {e}\n")
                        attempts += 1



#sys.exit()

if FLAG_search_gen:
    #2.) make list of genera where some species are missing but some already exists
    genera_dict =dict()
    with open(file_missing_but_genus_new, "r") as f:
        lines = f.readlines()[1:]
        for line in lines:
            line = line.strip() # remove whtie space and newline character from beginning and the end
            line = line.replace("\"","")
            line_information = line.split(",")
            if line_information[0] in genera_dict:
                genera_dict[line_information[0]].append(line_information[1])  
            else:
                genera_dict[line_information[0]] = [line_information[1]]


    #3.) make list of genera to search for
    genera_search = []
    with open(file_still_missing, "r") as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip() # remove whtie space and newline character from beginning and the end
                line = line.replace("\"","")
                line_information = line.split(",")
                if line_information[0] not in genera_search:
                    genera_search.append(line_information[0] )

    # 4.) check genera on ncbi and download
    with open(file_missing_genera, "w") as f_incomp, open(file_missing_genera2, "w") as f_miss, open(file_LOG, "a") as f_log, open(file_asseccion_info, "a") as f_assec, open(file_LOG2, "w") as f_log2:
        for genus in genera_search:
            print("")
            print(genus)
            #create search link
            new_search=""
            if genus in genera_dict:
                new_search = genus + "[ORGANISM] NOT (" +  "[ORGANISM] OR ".join(genera_dict[genus]) + "[ORGANISM]) " + search
                #print(new_search)
                #sys.exit()
            else:
                new_search = genus+ "[ORGANISM] " +  search

            attempts = 0
            while attempts < 3:
                try:
                    # search if genus is represented at ncbi except for already downloaded species
                    handle = Entrez.esearch(db=ncbi_db, retmax=50,term= new_search)
                    record = Entrez.read(handle)

                    #check if feedback is positive ore not
                    if len(record['IdList']) == 0:
                        print("Genus unknown: \t"+genus)
                        if genus in genera_dict:
                            f_miss.write(genus+","+ "&".join(genera_dict[genus])+ new_search +"\n") 
                        else:
                            f_incomp.write(genus+","+ new_search +"\n") 
                    else:
                        f_log2.write(genus+","+ new_search +"\n") 
                        #download genbank entries to get further organism info
                        handle = Entrez.efetch(db=ncbi_db, id=list(record['IdList']), rettype="genbank", retmode="text")
                        out_file=NCBI_folder + genus +"__NCBI.fasta"

                        #download fasta entries to save fasta file (otherwise empty files will be created)
                        handle2 = Entrez.efetch(db=ncbi_db, id=list(record['IdList']), rettype="fasta", retmode="text")
                        all_records = SeqIO.parse(handle2, 'fasta')
                        count = SeqIO.write(all_records, out_file, "fasta-2line")

                        # check if data are really created & 
                        # loop through all asseccion numbers and write into searchterm , organism, and file accession number 
                        if os.path.isfile(out_file):
                            for records in SeqIO.parse(handle, 'genbank'):
                                f_assec.write(genus+ "\t"+  records.annotations['organism']+"\t" + records.id + "\t" + out_file  +"\t" + (";").join(records.annotations['taxonomy']) + "\n")
       
                    break
                except:
                    print("Something went wrong")
                    attempts += 1
                    f_log.write(genus+","+ new_search +"\n") 

exit()

