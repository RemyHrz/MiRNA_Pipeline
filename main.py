#!/usr/bin/env python

import os, glob
import json
import math
import pandas as pd
from functools import reduce
import subprocess
from itertools import combinations
from Scripts.RNAnalysis_lib import write_tsv_from_dic,filter_table_with_ids,collapse_read,extract_targets,remove_common_target,keep_common_target,unique_readid_from_mirdp2,read_fasta,log_step,write_collapsed,write_collapsed_tab,same_directory,mirbase_uniq_org

# Setup working directory and paths

current_dir=os.getcwd()
main_dir=current_dir+"/RNA_Analysis"
script_path = os.path.dirname(os.path.realpath(__file__))+"/Scripts"

# creating necessary directories

dirs_list=["/cleanreads","/raw_fastq_files","/references/bowtie","/references/targets",
"/references/genomes","/references/mirbase","/miRDP2/collapsed_reads",
"/miRDP2/outputs","/results","/logs","/tables","/targeting"]

if os.path.exists(f"{main_dir}/config_default.json"):
    with open(f"{main_dir}/config_default.json", "r") as file:
        config=json.load(file)
    config["reads_dir"]=main_dir+"/cleanreads"
    config["raw_fastq_dir"]=main_dir+"/raw_fastq_files"
    config["logs_dir"]=main_dir+"/logs"
    config["bowtie_dir"]=main_dir+"/references/bowtie"
    config["genome_dir"]=main_dir+"/references/genomes"
    config["mirbase_dir"]=main_dir+"/references/mirbase"
    config["target_dir"]=main_dir+"/references/targets"
    config["results_dir"]=main_dir+"/results"
    config["miRDP2_dir"]=main_dir+"/miRDP2"
    config["tables_dir"]=main_dir+"/tables"
    config["targeting_dir"]=main_dir+"/targeting"
    config["collapsed_dir"]=config["miRDP2_dir"]+"/collapsed_reads"
    config["miRDP2_out_dir"]=config["miRDP2_dir"]+"/outputs"
    config["merged_genome_file"]=config["genome_dir"]+"/mergedGenome.fasta"
else:
    try:
        if not os.path.exists(main_dir):
            os.makedirs(main_dir)
        for directory in dirs_list:
            if not os.path.exists(main_dir+directory):
                os.makedirs(main_dir+directory)
        print("directories created, please add your files to the corresponding folder (raw_fastq_files,references/genomes,references/mirbase)")
# creating and loading config file
        config = {
    "quality_threshold":30,
    "mininum_rna_size":18,
    "maximum_rna_size":26,
    "trim_nt_front":4,
    "trim_nt_tail":4,
    "Adapter_sequence":"TGGAATTCTCGGGTGCCAAGGAACTCCAG",
    "threads":14,
    "bowtie_index_name":"reference",
    "collapsed_read_count_treshold":3,
    "mirbase_org_code":"ath",
    "fold_change_threshold":1,
    "dea_pvalue":0.05,
    "seed_bowtie":12345}
        with open(f"{main_dir}/config_default.json", "w") as file:
            json.dump(config, file)
        print("Default config file created, adjust the value as needed and restart the script")
        quit()
    except:
        print("directories and/or config file not created, there was a problem")
        quit()


#"""

main_log=f"{config["logs_dir"]}/logs"

with open(main_log, "w") as file:
    file.write("Starting Analysis")


# cleaning reads with fastp

log_step(f"{config["logs_dir"]}/logs","Cleaning reads ...")

cleanreads=f"{script_path}/cleanreads.sh {config["raw_fastq_dir"]} {config["reads_dir"]} {config["quality_threshold"]} {config["mininum_rna_size"]} {config["maximum_rna_size"]} {config["trim_nt_front"]} {config["trim_nt_tail"]} {config["logs_dir"]} {config["Adapter_sequence"]} {math.floor(config["threads"]/2)} 2>>{main_log}"
subprocess.run(cleanreads, shell = True, executable="/bin/bash")

log_step(f"{config["logs_dir"]}/logs","Done")

# creating a fasta file containing all reference genomes and indexing the resulting fasta file with bowtie1

if not os.path.exists(f"{config["bowtie_dir"]}/{config["bowtie_index_name"]}.1.ebwt"):
    log_step(f"{config["logs_dir"]}/logs","bowtie index not found, creating ...")
    with open(f"{config["genome_dir"]}/mergedGenome.fasta","a") as merged:
        for filename in os.listdir(config["genome_dir"]):
            if filename != "mergedGenome.fasta":
                with open(f"{config["genome_dir"]}/{filename}","r") as file:
                    merged.write(file.read())
    bowtie_build=f"bowtie-build {config["merged_genome_file"]} {config["bowtie_dir"]}/{config["bowtie_index_name"]} --threads {config["threads"]} 2>>{main_log} 1>>{main_log}"
    subprocess.run(bowtie_build, shell = True, executable="/bin/bash")
    log_step(f"{config["logs_dir"]}/logs","Done")

log_step(f"{config["logs_dir"]}/logs","Bowtie mapping ...")

# Mapping each cleaned fastq file to the indexed mergedgenome with bowtie

for filename in os.listdir(config["reads_dir"]):
    outname=os.path.basename(filename).split(".")[0]+f"on{config["bowtie_index_name"]}.sam"
    log_step(f"{config["logs_dir"]}/logs",outname)
    bowtie_mapping=f"bowtie --seed {config["seed_bowtie"]} -p {config["threads"]} -x {config["bowtie_dir"]}/{config["bowtie_index_name"]} -q {config["reads_dir"]}/{filename} -S {config["results_dir"]}/{outname} 2>> {main_log}"
    subprocess.run(bowtie_mapping, shell = True, executable="/bin/bash")
    
log_step(f"{config["logs_dir"]}/logs","Done")

# Collapsing read : removing redundacy while keeping the occurence of each sequence, creating tables for each 
# fastq file containing collapsed read_id, Sequence and occurence of said sequence in the fastq file

log_step(f"{config["logs_dir"]}/logs","collapsing reads for miRDP2 ...")

for file in os.listdir(config["reads_dir"]):
    basename=os.path.basename(file).split(".")[0]
    outname=basename+f".collapsed.x{config["collapsed_read_count_treshold"]}.fa"
    log_step(f"{config["logs_dir"]}/logs",outname)
    read_dict=collapse_read(f"{config["reads_dir"]}/{file}", config["collapsed_read_count_treshold"])
    write_collapsed(read_dict,f"{config["collapsed_dir"]}/{outname}")
    write_collapsed_tab(read_dict,f"{config["tables_dir"]}/{basename}.tab",basename)

log_step(f"{config["logs_dir"]}/logs","Done")

# Identification of known and novel miRNAs with miRDP2

log_step(f"{config["logs_dir"]}/logs","Identification of known and potential novel miRNAs with miRDP2 ...")

for filename in os.listdir(config["collapsed_dir"]):
    log_step(f"{config["logs_dir"]}/logs",f"{os.path.basename(filename).split(".fa")[0]}on{config["bowtie_index_name"]}")
    miRDP2=f"miRDP2-v1.1.4_pipeline.bash -g {config["merged_genome_file"]} -x {config["bowtie_dir"]}/reference -f -i {config["collapsed_dir"]}/{filename} -o {config["miRDP2_out_dir"]} -p {config["threads"]} 2>>{config["logs_dir"]}/logs"
    subprocess.run(miRDP2, shell = True, executable="/bin/bash")

log_step(f"{config["logs_dir"]}/logs","Done")

log_step(f"{config["logs_dir"]}/logs","Extracting miRDP2 results ...")

# Creating lists of resulting files of interest from all miRDP2 subdirectories

prediction_files_list=[]
for file in glob.glob(f"{config["miRDP2_out_dir"]}/*/*_filter_P_prediction"):
    prediction_files_list.append(file)

known_files_list=[]
for file in glob.glob(f"{config["miRDP2_out_dir"]}/*/known_miR.aln"):
    known_files_list.append(file)

# Matching pair of file from the same directory

matched_pairs = [(a, b) for a, b in zip(known_files_list, prediction_files_list) if same_directory(a, b)]

log_step(f"{config["logs_dir"]}/logs","Merging results...")

# Extracting read ids from the results of miRDP2 and retrieving corresponding 
# sequence from tables of collapsed read 

df_list=[]
for pair in matched_pairs:
    read_ids=unique_readid_from_mirdp2(pair[0],pair[1])
    for table in os.listdir(config["tables_dir"]):
        exp_name=os.path.basename(table).split(".")[0]
        if exp_name in str(pair[0]):
            df=filter_table_with_ids(f"{config["tables_dir"]}/{table}",read_ids)
            break
    df_list.append(df)

# merging sequences from each replicate into a one-column dataframe to gather 
# all sequences while removing redundancy

merged_df = reduce(lambda left, right: pd.merge(left, right, on=["Sequence"], how="outer").fillna(0), df_list)
merged_df.to_csv(f"{config["results_dir"]}/seqtable.tmp", sep="\t", index=False)

# merging occurrence table of each replicates to the sequence df to creating 
# an occurence table of each sequence in each replicate

for file in os.listdir(config["tables_dir"]):
    df=pd.read_csv(f"{config["tables_dir"]}/{file}", sep="\t")
    df=df.drop("read_id",axis=1)
    merged_df=pd.merge(merged_df, df, on=["Sequence"], how="left").fillna(0)

merged_df.to_csv(f"{config["results_dir"]}/miRtable.tmp", sep="\t", index=False)

# writing list of all sequences into a fasta file to reidentify miRNA using bowtie, as miRDP2 
# can assign the same read sequence to different miRNA that share the same sequence

log_step(f"{config["logs_dir"]}/logs","Done")

Sequences=merged_df["Sequence"]

with open(f"{main_dir}/all_seq.fa", "w") as fasta:
    count_id=0
    for seq in Sequences:
        fasta.write(f">seq{count_id}\n{seq}\n")
        count_id+=1

log_step(f"{config["logs_dir"]}/logs","Mapping of all sequence against reference for identification of miRNAs")

# deduplicating sequence from mirbase reference and selecting mirna of the selected organism

mirbase_dic=mirbase_uniq_org(f"{config["mirbase_dir"]}/mature.fa",config["mirbase_org_code"])

with open(f"{config["mirbase_dir"]}/mirbase_org_uniq.fa","w") as file:
    for seq,mirna_id in mirbase_dic.items():
        file.write(f">{mirna_id}\n{seq}\n")

# building mirbase index for bowtie

bowtie_build_mirbase=f"bowtie-build {config["mirbase_dir"]}/mirbase_org_uniq.fa {config["bowtie_dir"]}/mirbase_uniq --threads {config["threads"]} 2>>{main_log} 1>>{main_log}"

subprocess.run(bowtie_build_mirbase, shell = True, executable="/bin/bash")
    
# mapping all mirna sequence for identification

bowtie_mirbase=f"bowtie --seed {config["seed_bowtie"]} --best --strata -k 1 -m 1 -p {config["threads"]} -x {config["bowtie_dir"]}/mirbase_uniq -f {main_dir}/all_seq.fa 1>{config["results_dir"]}/knownmir.aln 2>> {main_log}"

subprocess.run(bowtie_mirbase, shell = True, executable="/bin/bash")

log_step(f"{config["logs_dir"]}/logs","Done")

log_step(f"{config["logs_dir"]}/logs","Merging results of mapping and miRDP2 ...")

# merging known mirna names from bowtie result to the occurence table,

aln=pd.read_csv(f"{config["results_dir"]}/knownmir.aln", sep="\t", header=None,usecols=[0,2])

aln.columns = ["Sequence", "miR_Name"]

aln["Sequence"]=aln["Sequence"].map(read_fasta(f"{main_dir}/all_seq.fa"))

Table=pd.merge(merged_df, aln, on=["Sequence"], how="left").fillna("unknown")
column_order = ["Sequence","miR_Name"] + [col for col in Table.columns if col not in ["Sequence","miR_Name"]]
Table = Table.loc[:, column_order]

# ensuring occurence per replicate columns are integer

for col in Table.columns[2:]:
    Table[col] = Table[col].astype(int)

# adding  a rolling number to unknwon mirna in order to discriminate them

counter = iter(range(1, len(Table[Table['miR_Name'] == 'unknown']) + 1))
Table["miR_Name"]=Table["miR_Name"].apply(lambda x: f"unknown_miR{next(counter)}" if x == "unknown" else x)

# ensuring mir_name and Sequence are strings

for col in Table.columns[0:1]:
    Table[col]= Table[col].astype(str)

Table.to_csv(f"{config["results_dir"]}/miRtable.tsv", sep="\t", index=False)

log_step(f"{config["logs_dir"]}/logs","Done")

# performing differential expression analysis with DEseq2 in rscript 

log_step(f"{config["logs_dir"]}/logs","Differential expression analysis with DESeq2")

subprocess.run(["Rscript", "--vanilla", "--quiet",f"{script_path}/deseq.R",f"{main_dir}",f"{config["fold_change_threshold"]}",f"{config["dea_pvalue"]}"], check=True)

log_step(f"{config["logs_dir"]}/logs","Done")


# Target analysis of each replicate Differentially expressed mirna on 
# each target file with target finder 

log_step(f"{config["logs_dir"]}/logs","Targeting ...")

for file in glob.glob(f"{config["targeting_dir"]}/*.tsv"):
    exp_name=os.path.basename(file).split('.')[0]+"_"
    for target in glob.glob(f"{config["target_dir"]}/*.fasta"):
        outname=exp_name+os.path.basename(target).split('.')[0]
        log_step(f"{config["logs_dir"]}/logs",f"{outname}")
        targetingshell=f"{script_path}/targeting.sh {file} {target} {math.floor((config["threads"]-1)/2)} >>{config["targeting_dir"]}/{outname}.target"
        subprocess.run(targetingshell, shell = True, executable="/bin/bash")

log_step(f"{config["logs_dir"]}/logs","Done")

# Extracting targets by experimental condition

log_step(f"{config["logs_dir"]}/logs","Extracting target list ...")

#"""

targeted_genes_dic={}
for file in glob.glob(f"{config["targeting_dir"]}/*.target"):
    filename=os.path.basename(file)
    name_parts=filename.split("_")
    virus=name_parts[0]
    regulation=name_parts[1]
    targetfile="_".join(name_parts[2:-1])+name_parts[-1].split(".")[0]
    keyname=f"{virus}-{regulation}-{targetfile}"
    targeted_genes_dic[keyname]=file
    
targeted_genes_dic_processed={}
for key1,key2 in combinations(targeted_genes_dic.keys(),2):
    keys1=key1.split("-")
    keys2=key2.split("-")
    # compare down and upregulated of same virus & same file
    if (keys1[0]==keys2[0]) & (keys1[-1]==keys2[-1]) & (keys1[1] != keys2[1]):
        target1=extract_targets(targeted_genes_dic[key1])
        target2=extract_targets(targeted_genes_dic[key2])
        keyname_up=f"{keys1[0]}_UP_{keys1[-1]}"
        keyname_down=f"{keys1[0]}_DOWN_{keys1[-1]}"
        if "-UP-" in key1:
            targeted_genes_dic_processed[keyname_up]=remove_common_target(target1,target2)
            targeted_genes_dic_processed[keyname_down]=remove_common_target(target2,target1)
        if "-DOWN-" in key1:
            targeted_genes_dic_processed[keyname_down]=remove_common_target(target1,target2)
            targeted_genes_dic_processed[keyname_up]=remove_common_target(target2,target1)
    # compare same regulation of different virus
    if (keys1[0]!=keys2[0]) & (keys1[-1]==keys2[-1]) & (keys1[1] == keys2[1]):
        target1=extract_targets(targeted_genes_dic[key1])
        target2=extract_targets(targeted_genes_dic[key2])
        if keys1[1]=="UP":
            keyname=f"{keys1[0]}_{keys2[0]}_UP_{keys1[-1]}"
        else:
            keyname=f"{keys1[0]}_{keys2[0]}_DOWN_{keys1[-1]}"
        targeted_genes_dic_processed[keyname]=keep_common_target(target1,target2)

write_tsv_from_dic(read_fasta(f"{config["mirbase_dir"]}/mirbase_org_uniq.fa"),f"{config["targeting_dir"]}/mirbase.tab")
targetingshell=f"{script_path}/targeting.sh {file} {target} {math.floor((config["threads"]-1)/2)} >{config["targeting_dir"]}/{outname}.target"
subprocess.run(targetingshell, shell = True, executable="/bin/bash")
 
df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in targeted_genes_dic_processed.items()]))
df = df.fillna("")

df.to_csv(f"{config["results_dir"]}/targets.tsv",sep="\t", index=False)

log_step(f"{config["logs_dir"]}/logs","Done")

log_step(f"{config["logs_dir"]}/logs","Removing temp file from result directory")

if os.path.exists(f"{config["results_dir"]}/seqtable.tmp"):
    os.remove(f"{config["results_dir"]}/seqtable.tmp")

if os.path.exists(f"{config["results_dir"]}/miRtable.tmp"):
    os.remove(f"{config["results_dir"]}/miRtable.tmp")

log_step(f"{config["logs_dir"]}/logs","Analysis finished")


