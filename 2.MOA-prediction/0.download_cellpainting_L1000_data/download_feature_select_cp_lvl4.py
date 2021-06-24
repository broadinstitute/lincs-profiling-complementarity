import os
import argparse
import pathlib
import pandas as pd
import numpy as np
from collections import defaultdict
from pycytominer import feature_select


def recode_dose(dose_value):
    """This function recode the doses in Level-4 data to 8 distinct dose classes"""
    
    doses = [0.04,0.12,0.37,1.11,3.33,10.0,20.0,25.0]
    for x in range(len(doses)-1):
        if (dose_value > 0.0) & (dose_value <= 0.04):
            dose_value = 0.04
        elif doses[x] <= round(dose_value,2) < doses[x+1]:
            dose_value = doses[x]
    return dose_value
    
def feature_selection(df_lvl4): 
    """
    Perform feature selection by dropping columns with null values 
    (greater than 384 i.e. equivalent to one plate worth of cell profiles) 
    and highly correlated values from the data.
    """
    metadata_columns = [x for x in df_lvl4.columns if (x.startswith("Metadata_"))]
    df_lvl4_metadata = df_lvl4[metadata_columns].copy()
    df_lvl4_features = df_lvl4.drop(metadata_columns, axis = 1)
    null_cols = [col for col in df_lvl4_features.columns if df_lvl4_features[col].isnull().sum() > 384]
    df_lvl4_features.drop(null_cols, axis = 1, inplace=True)
    ##feature selection was done already..prior to getting the spherized data!!
    ###df_lvl4_features = feature_select(df_lvl4_features, operation=["correlation_threshold", "variance_threshold"])
    
    for col in df_lvl4_features.columns:
        if df_lvl4_features[col].isnull().sum():
            df_lvl4_features[col].fillna(value=df_lvl4_features[col].mean(), inplace = True)
            
    df_meta_info = df_lvl4_metadata[['Metadata_broad_sample', 'Metadata_pert_id', 'Metadata_Plate', 'Metadata_Well',
                                     'Metadata_broad_id', 'Metadata_moa', 'Metadata_dose_recode']].copy()
    df_lvl4_new = pd.concat([df_meta_info, df_lvl4_features], axis=1)
    
    return df_lvl4_new
    
def merge_dataframe(df, pertinfo_file):
    """
    This function merge aligned L1000 and Cell painting Metadata information dataframe 
    with the Level-4 data, change the values of the Metadata_dose_recode column 
    and create a new column 'replicate_name' that represents each replicate in the dataset
    """ 
    df_pertinfo = pd.read_csv(pertinfo_file)
    df_lvl4_new = df.merge(df_pertinfo, on='Metadata_broad_sample', how = 'outer')
    no_cpds_df = df_lvl4_new[df_lvl4_new['pert_iname'].isnull()].copy().reset_index(drop = True)
    df_lvl4_new.drop(df_lvl4_new[df_lvl4_new['pert_iname'].isnull()].index, inplace = True)
    df_lvl4_new.reset_index(drop= True, inplace = True)
    df_lvl4_new['Metadata_dose_recode'] = df_lvl4_new['Metadata_dose_recode'].map({0.0:0,0.04:1,0.12:2,0.37:3,1.11:4,
                                                                                   3.33:5,10.0:6,20.0:7})
    df_lvl4_new['replicate_name'] = ['replicate_' + str(x) for x in range(df_lvl4_new.shape[0])]
    
    return df_lvl4_new, no_cpds_df
    
def save_to_csv(df, path, file_name, compress=None):
    """saves dataframes to csv"""
    
    if not os.path.exists(path):
        os.mkdir(path)
    
    df.to_csv(os.path.join(path, file_name), index=False, compression=compress)
    
def download_cp_lvl4_align_moa(profile_link=None, pertinfo_file=None, data_dir=None):
    """
    This function downloads the Cell painting level-4 profiles data from its github repo, 
    perform dose recode, feature selection and align its compounds and MOAs (Mechanism of actions) 
    with the ones found in L1000 profiles and save the resulting dataframe as csv file.
    Args:
            profile_link: Github repo link of Cell painting level-4 data.
            pertinfo_file: Aligned metadata perturbation csv file directory for Cell painting & L1000
            data_dir: Directory to save the cell painting level-4 data.
            
    Output:
            saved csv file: Cell painting level-4 profiles data
    """
    data_dir = 'D:\cell_painting_profiles\profiles\cellpainting_lvl4_cpd_replicate_datasets'
    df_level4 = pd.read_csv(profile_link, compression='gzip',low_memory = False)
    df_level4['Metadata_dose_recode'] = df_level4['Metadata_mmoles_per_liter'].apply(recode_dose)
    df_level4_new = feature_selection(df_level4)
    df_level4_new, df_level4_no_cpds = merge_dataframe(df_level4_new, pertinfo_file)
    save_to_csv(df_level4_new, data_dir, 'cp_level4_cpd_replicates.csv.gz', compress="gzip")
    print("Done!! downloaded cell-painting level-4 profiles data and cleaned it for further analysis")

def parse_args():
    """Arguments to pass to this Script"""
    
    parser = argparse.ArgumentParser(description="Parse arguments")
    parser.add_argument('profile_link', type=str, default = "https://github.com/broadinstitute/lincs-cell-painting/blob/e17a47c9a524d4789982511dd5db9b0202ff6cc8\
/spherized_profiles/profiles/2016_04_01_a549_48hr_batch1_dmso_spherized_profiles_with_input_normalized_by_whole_plate.csv.gz?raw=true",
                        nargs='?', help='Github repo link of Cell painting level-4 data')
    parser.add_argument('pertinfo_file', type=str, default = '../aligned_moa_CP_L1000.csv', nargs='?', help='Aligned metadata perturbation.csv file directory for Cell painting & L1000')
    parser.add_argument('data_dir', type=str, default = os.getcwd(), nargs='?', help='Directory to save the cell painting \
    level-4 data')
    return parser.parse_args()
    
if __name__ == '__main__':
    args = parse_args()
    download_cp_lvl4_align_moa(args.profile_link, args.pertinfo_file, args.data_dir)