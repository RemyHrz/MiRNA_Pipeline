#!/usr/bin/env python

import os
import pandas as pd
import random
import datetime

def read_fasta(fasta)-> dict:
    fasta_dict = {}
    with open(fasta) as f:
        name = None
        for line in f:
            line = line.strip()
            if line.startswith(">"):  # Header line
                name = line[1:]
                fasta_dict[name] = ""
            else:
                fasta_dict[name] = line
    return fasta_dict

def log_step(main_log,step:str):
    with open(main_log, "a") as file:
        file.write(f"\n{datetime.datetime.now().strftime('%d/%m/%y %H:%M:%S')} {step}\n")
    print(f"{step}\n")

def collapse_read(input_file,threshold:int):
    read_dict={}
    with open(input_file,"r") as file:
        line=file.readline()
        while line:
            line = file.readline().strip()
            if line not in read_dict:
                read_dict[line] = 1
            else:
                read_dict[line] +=1
            line=file.readline()
            line=file.readline()
            line=file.readline()
    clean_dict={}
    for seq,count in read_dict.items():
            if count >= threshold:
                clean_dict[seq]=count
    return clean_dict

def write_collapsed(read_dict:dict,output_file):
    read_num=1
    with open(output_file,"w") as file:
        for seq,count in read_dict.items():
                file.write(f">read{read_num}_x{count}\n{seq}\n")
                read_num+=1

def same_directory(path1, path2):
    return os.path.dirname(path1) == os.path.dirname(path2)

def unique_readid_from_mirdp2(known_file,predict_file):
    col_a = pd.read_csv(known_file,sep="\t", usecols=[0], header=None).squeeze("columns")
    col_b = pd.read_csv(predict_file,sep="\t", usecols=[2], header=None).squeeze("columns")
    combined = pd.concat([col_a, col_b], ignore_index=True)
    unique_readid = combined.unique()
    return unique_readid

def write_collapsed_tab(read_dict:dict,output_file,exp_name):
    with open(output_file,"w") as file:
        read_num=1
        file.write(f"Sequence\t{exp_name}\tread_id\n")
        for seq,count in read_dict.items():
                file.write(f"{seq}\t{count}\tread{read_num}_x{count}\n")
                read_num+=1

def filter_table_with_ids(table,read_ids):
    df=pd.read_csv(table, sep="\t")
    filtered_df = df[df['read_id'].isin(read_ids)]
    return filtered_df["Sequence"]

def mirbase_uniq_org(input_file,mirbase_org_code):
    fasta_dict = {}
    with open(input_file) as f:
        while True:
            mirna_id = f.readline().strip().split(" ")[0]
            seq = f.readline().strip().replace("U","T")
            if not mirna_id:
                break
            if mirna_id.startswith(">") and mirbase_org_code in mirna_id: # Header line
                name = mirna_id[1:]
                fasta_dict[seq] = name
    return fasta_dict

def extract_targets(file):
    df=pd.read_csv(file,header=None,sep="|",usecols=[0])
    targets_ids=df[0].str.split("\t",expand=True)
    targets_ids=targets_ids[1].str.replace(r"\..*", "", regex=True)
    return targets_ids.unique()

def remove_common_target(lista,listb):
    seta=set(lista)
    setb=set(listb)
    return list(seta - setb)

def keep_common_target(lista,listb):
    seta=set(lista)
    setb=set(listb)
    return list(seta & setb)

def write_tsv_from_dic(read_dict:dict,output_file):
    with open(output_file,"w") as file:
        for seq_id,seq in read_dict.items():
                file.write(f"{seq_id}\t{seq}\n")




