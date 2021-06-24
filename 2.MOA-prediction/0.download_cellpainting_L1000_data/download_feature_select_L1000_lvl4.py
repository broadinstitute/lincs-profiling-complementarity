import os
import requests
import pickle
import argparse
import pandas as pd
import numpy as np
import re
from os import walk
from collections import Counter
import cmapPy.pandasGEXpress.parse_gct as pg
from cmapPy.pandasGEXpress.parse import parse
from io import BytesIO
from urllib.request import urlopen
from zipfile import ZipFile

def construct_lvl4_data(data_dir, level4_dir, pertinfo_file):
    """
    This function returns L1000 Level-4 Data that is aligned with 
    the important metadata information to compute compound's median scores
    """
    
    lvl4_data = parse(os.path.join(data_dir, level4_dir))
    lvl4_data = lvl4_data.data_df.rename_axis(None).T
    lvl4_data = lvl4_data.rename_axis(None).reset_index()
    lvl4_data.rename(columns={"index": "replicate_id"}, inplace = True)
    df_metalvl_5 = pd.read_csv(os.path.join(data_dir, 'col_meta_level_5_REP.A_A549_only_n9482.txt'), delimiter = "\t")
    lvl4_data['sig_id'] = lvl4_data['replicate_id'].apply(lambda x: ''.join(re.sub('(24H.+(\_|\.)[A-Z0-9]+.)\:', '24H:', x)))
    df_meta_features = df_metalvl_5[['sig_id', 'pert_id', 'pert_idose', 'det_plate', 'det_well']].copy()
    df_meta_features['dose'] = df_meta_features['pert_idose'].map({'-666' : 0, '0.04 uM' : 1, '0.12 uM' : 2, '0.37 uM' : 3,
                                                                   '1.11 uM' : 4, '3.33 uM' : 5, '10 uM' : 6, '20 uM' : 7})
    df_pertinfo = pd.read_csv(pertinfo_file)
    df_pertinfo.rename(columns={"broad_id": "pert_id"}, inplace = True)
    df_meta_features = pd.merge(df_meta_features, df_pertinfo, on=['pert_id'], how = 'left')
    lvl4_data = pd.merge(lvl4_data, df_meta_features, on='sig_id')
    
    return lvl4_data
    
def feature_selection(df_data):
    
    """
    Perform feature selection by dropping columns with null compounds values, 
    and highly correlated landmark genes from the data.
    """
    
    df_data_genes = df_data.drop(['replicate_id', 'Metadata_broad_sample', 'pert_id', 'dose', 'pert_idose', 
                                  'pert_iname', 'moa', 'sig_id', 'det_plate', 'det_well'], axis = 1).copy()
    df_data_corr = df_data_genes.corr(method = 'spearman')
    drop_cols = []
    n_cols = len(df_data_corr.columns)
    for i in range(n_cols):
        for k in range(i+1, n_cols):
            val = df_data_corr.iloc[k, i]
            col = df_data_corr.columns[i]
            if abs(val) >= 0.9:
                drop_cols.append(col)
    df_data.drop(set(drop_cols), axis = 1, inplace = True)
    df_data.drop(df_data[df_data['pert_iname'].isnull()].index).reset_index(drop = True, inplace = True)
    
    return df_data

def save_to_csv(df, path, file_name, compress=None):
    """saves dataframes to csv"""
    
    if not os.path.exists(path):
        os.mkdir(path)
    
    df.to_csv(os.path.join(path, file_name), index=False, compression=compress)

def download_L1000_lvl4_align_moa(L1000_datalink=None, pertinfo_file=None, data_dir=None):
    """
    This function downloads the L1000 data from figshare, perform feature selection on L1000 Level-4 profiles 
    and align its compounds and MOAs (Mechanism of actions) with the ones found in Cell painting profiles and 
    save the resulting dataframe as csv file.
    Args:
            profile_link: figshare link of L1000 data.
            pertinfo_file: Aligned metadata perturbation csv file directory for Cell painting & L1000
            data_dir: Directory where all the L1000 data are saved.
            
    Output:
            Directory: L1000 data
    """
    zipurl = L1000_datalink
    ##download from figshare
    if not os.path.exists(data_dir):
        os.mkdir(data_dir)
    with urlopen(zipurl) as zipresp:
        with ZipFile(BytesIO(zipresp.read())) as zfile:
            zfile.extractall(data_dir)
            
    df_level4 = construct_lvl4_data(data_dir, 'level_4_zspc_n27837x978.gctx', pertinfo_file)
    df_level4 = feature_selection(df_level4)
    save_to_csv(df_level4, data_dir, 'L1000_level4_cpd_replicates.csv.gz', compress="gzip")
    print("Done!! downloaded L1000 data and cleaned the Level-4 profiles data for further analysis")
    
def parse_args():
    """Arguments to pass to this Script"""
    
    parser = argparse.ArgumentParser(description="Parse arguments")
    parser.add_argument('--L1000_datalink', type=str, default = "https://ndownloader.figshare.com/articles/13181966/versions/1", nargs='?', 
                        help='L1000 data repository link')
    parser.add_argument('--pertinfo_file', type=str, default = '../aligned_moa_CP_L1000.csv', nargs='?', 
                        help='Aligned metadata perturbation csv file directory for Cell painting & L1000')
    parser.add_argument('--data_dir', type=str, default = os.getcwd(), nargs='?', 
                        help='Directory where all the L1000 data are saved.')
    return parser.parse_args()
    
if __name__ == '__main__':
    args = parse_args()
    download_L1000_lvl4_align_moa(args.L1000_datalink, args.pertinfo_file, args.data_dir)