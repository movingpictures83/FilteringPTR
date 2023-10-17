#!/usr/bin/env python
# coding: utf-8

# # Supplementary Notebook 1: Filtering PTR values
# ## Paper: Novel Approach for Microbiome Analysis Using Bacterial Replication Rates and Causal Inference to Determine Resistome Potential
# ### Vitalii Stebliankin, Musfiqur Sazal, Camilo Valdes, Kalai Mathee, and GiriNarasimhan
# 
# #### Dataset: Gibson et al. (BioProject ID: PRJNA301903)
# 
# ### 1. Explore PTR Distribution

# In[1]:


#create output folder
import os

import PyPluMA
import PyIO
# Read the PTR data:
import pandas as pd
# Visualize the data distribution
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def merge_strains(ptr_df, taxonomy_file, abundance_path, out_path):
    ################################
    # Merge strains to species level
    ################################

    samples_list = list(ptr_df.columns)
    samples_list.remove("NZ")
    samples_list = [x.strip("\n") for x in samples_list]

    ptr_df.rename(columns=lambda x: x.strip("\n") + "_PTR", inplace=True)
    ptr_df.rename(columns={"NZ_PTR": "NZ"}, inplace=True)
    # Step 1 - Merge with abundance file
    abundance_df = pd.read_csv(abundance_path)

    abundance_df.rename(columns=lambda x: x.strip("\n") + "_abundance", inplace=True)
    abundance_df.rename(columns={"NZ_abundance": "NZ"}, inplace=True)

    nz_list = list(ptr_df["NZ"])
    abundance_df = abundance_df[abundance_df["NZ"].isin(nz_list)]
    # ptr_df = ptr_df.fillna(0) # convenient for computation. Note: Change for 1 later
    abundance_df = abundance_df.fillna(0)
    abundance_df = abundance_df.reset_index(drop=True)

    ptr_df = ptr_df.merge(abundance_df, on="NZ", how="left")

    # Step 2 - Merge NZ with taxa species name
    taxonomy_df = pd.read_csv(taxonomy_file, sep="\t")

    taxonomy_df["NZ"] = taxonomy_df["NZ"].apply(lambda x: x.split(".")[0])
    taxonomy_cols = list(taxonomy_df.columns)
    ptr_df = taxonomy_df.merge(ptr_df, on="NZ", how="right")

    species_list = list(ptr_df["species"].unique())

    # Step 3 - Combine taxa

    # create dictionary with new values
    ptr_dict = {}
    ptr_dict["sample"] = samples_list

    previous_species = ""

    k = -1
    columns_list = []
    for specie in species_list:
        print("Merging " + specie + "...")
        if isinstance(specie, str):
            tmp_df = ptr_df[ptr_df["species"] == specie]
            ptr_dict[specie + "#abundance"] = []
            ptr_dict[specie + "#PTR"] = []
            for sample in ptr_dict["sample"]:
                tmp_df_sample = tmp_df[["NZ", "species", sample + "_PTR", sample + "_abundance"]]
                if len(tmp_df_sample)>0:
                    total_abundance_not_droppped = tmp_df_sample[sample + "_abundance"].sum()
                else:
                    total_abundance_not_droppped=0
                tmp_df_sample = tmp_df_sample.dropna()
                total_abundance = tmp_df_sample[sample + "_abundance"].sum()

                if len(tmp_df_sample) > 0:
                    tmp_df_sample[sample + "_cumulative"] = tmp_df.apply(
                        lambda row: row[sample + "_abundance"] * row[sample + "_PTR"] / total_abundance, axis=1)
                    ptr = tmp_df_sample[sample + "_cumulative"].sum()
                    ptr_dict[specie + "#abundance"].append(total_abundance_not_droppped)
                    ptr_dict[specie + "#PTR"].append(ptr)

                else:
                    ptr_dict[specie + "#abundance"].append(total_abundance_not_droppped)
                    ptr_dict[specie + "#PTR"].append(np.nan)
        columns_list.append(specie + "#abundance")
        columns_list.append(specie + "#PTR")

    new_ptr_df = pd.DataFrame.from_dict(ptr_dict)
    new_ptr_df = new_ptr_df[["sample"]+ columns_list]
    new_ptr_df.to_csv(out_path, index=False)
    return new_ptr_df

def filter_major(cutoff, out_file, ptr_df, species_list, metadata_file):
    # Filter out the species that present in less than "cutoff" samples:
    updated_n_sp=0
    for sp in species_list:
        tmp_df = ptr_df[ptr_df[sp+"#PTR"].notnull()]
        if len(tmp_df)<cutoff:
            ptr_df = ptr_df.drop(sp+"#PTR", axis=1)
            ptr_df = ptr_df.drop(sp + "#abundance", axis=1)
        else:
            updated_n_sp+=1
        # sns.regplot(data=tmp_df, x=sp+"#PTR", y=sp+"#abundance")
        # plt.show()
    pass
    pass
    print("Number of selected species: {}".format(updated_n_sp))
    # Clinical variables:
    metadata_df = pd.read_csv(metadata_file)

    merged_df = metadata_df.merge(ptr_df, on="sample", how="right")
    merged_df = merged_df[merged_df["Day_of_Life"].notnull()]
    merged_df = merged_df.drop_duplicates()
    merged_df.to_csv(out_file, index=False)

#    merged_df = merged_df.fillna(1)
    merged_df.to_csv(out_file, index=False)
    return merged_df

def combine(row, species_list):
    total_abundance = 0
    cum_ptr = 0
    for species in species_list:
        #try:
        ptr = row[species+"#PTR"]
        abundance = row[species+"#abundance"]
        if (not np.isnan(ptr)) and (not np.isnan(abundance)):
            cum_ptr += ptr*abundance
            total_abundance += abundance
            pass
            pass
    if total_abundance==0:
        return np.nan
    else:
        return cum_ptr/total_abundance

class FilteringPTRPlugin:
 def input(self, inputfile):
    self.parameters = PyIO.readParameters(inputfile)
 def run(self):
    pass
 def output(self, outputfile):
    out_dir = outputfile
    #out_dir = "analysis-out/1-FilteringPTR"
    #intermediate_dir = "{}/intermediate_files".format(out_dir)
    #if not os.path.exists(out_dir):
    #os.mkdir(out_dir)
    #if not os.path.exists(intermediate_dir):
    #os.mkdir(intermediate_dir)

    # Reading PTR file for all bacteria
    ptr_df = pd.read_csv(PyPluMA.prefix()+"/"+self.parameters["mergedPTR"], index_col=0)
    #ptr_df = pd.read_csv("A-out/merged_ptr.csv", index_col=0)
    ptr_df = ptr_df.T
    ptr_df

    # Getting list of all genomes:
    all_genomes = list(ptr_df.columns)

    # Get the list of all PTR values:
    all_ptr_list = []
    for genome in all_genomes:
       all_ptr_list+= list(ptr_df[genome].dropna())

    # Plot the distribution:
    ax = sns.distplot(all_ptr_list)
    plt.title("PTR Distribution")


    # Some of the PTR values are as larg as value 14. However, it is unlikely that any microbial subpopulation has on average 14 replication forks. Therefore, 2% larges PTR values were considered as outliers.
    # 
    # ### 2 - Filter PTR outliers

    # In[2]:

    import numpy as np
    drop_persent = 2
    #Drop 2% largest PTR
    all_ptr_list.sort() # sort in increasing order
    top_ptr_values_index = int((len(all_ptr_list)/100) * (100-drop_persent))
    cutoff_value = all_ptr_list[top_ptr_values_index-1]
    all_ptr_list_prunned = all_ptr_list[0:top_ptr_values_index]

    # Remove 2% largest from PTR Data Frame:
    prunned_ptr_df = ptr_df.copy()
    prunned_ptr_df[all_genomes] = ptr_df[all_genomes].applymap(lambda x: x if x<cutoff_value else np.nan)

    # Save filtered dataframe:
    prunned_ptr_df.to_csv(PyPluMA.prefix()+"/"+self.parameters["intermediate"]+"/PTR_strain_filtered.csv")
    intermediate_dir = PyPluMA.prefix()+"/"+self.parameters["intermediate"]
    # Plot new distribution:
    ax = sns.distplot(all_ptr_list_prunned)
    plt.title("PTR distribution after dropping {}% largest PTR values".format(drop_persent))
    plt.show()


    # ### 3 - Averaging PTR to species level
    # Here we obtain the average number of replication forks per one genome in the species population.

    # In[3]:


    # Combine by species
    ptr_merged_path = intermediate_dir+"/PTR_strain_filtered.csv"
    taxonomy_file = PyPluMA.prefix()+"/"+self.parameters["taxonomyfile"]
    #taxonomy_file = "metadata/taxonomy.txt" # NCBI taxonomy. The file was obtained from RefSeqv.92. 
    abundance_path = PyPluMA.prefix()+"/"+self.parameters["mergedabundance"]
    #abundance_path = "A-out/merged_abundance.csv"

    out_path = intermediate_dir+"/PTR_species_filtered.csv"

    ptr_df = pd.read_csv(ptr_merged_path, index_col=0)
    ptr_df = ptr_df.T
    ptr_df["NZ"] = ptr_df.index
    ptr_df = ptr_df.reset_index(drop=True)
    ptr_df = merge_strains(ptr_df, taxonomy_file, abundance_path, out_path)


    # ### 4 - Compute Average PTR across samples

    # In[4]:


    ptr_file = intermediate_dir+"/PTR_species_filtered.csv"
    ptr_combined = intermediate_dir+"/PTR_sample_filtered.csv" #output
    metadata_path = PyPluMA.prefix()+"/"+self.parameters["metadatapath"]
    #metadata_path = "metadata/metadata_gibson.csv" #Metadata for Gibson et al. dataset

    df = pd.read_csv(ptr_file, index_col=0)

    relevant_coluns = PyIO.readSequential(PyPluMA.prefix()+"/"+self.parameters["relevantcolumns"])
    #relevant_coluns = ["sample", "Day_of_Life", "PostMenst_Age", "Individual", "Gestational_Age",
    #               "Birthweight", "Gentamicin", "Cefazolin","Ampicillin", "Trimethoprim-Sulfamathoxazole", "Meropenem",
    #               "Vancomycin", "Ticarcillin-Clavulanate", "Clindamycin", "Cefotaxime", "Total_abx", "r_Gentamicin",
    #               "r_Meropenem", "r_Ticarcillin-Clavulanate", "r_Vancomycin", "r_Ampicillin",
    #               "r_Cefotaxime","r_TOTAL","Human_Milk","Maternal_Milk", "Donor_Milk", "Formula","Fortification","Vitamin_A",
    #               "Caffeine","Iron","Furosemide_Lasix","m_ampicillin","m_ceftriaxone","m_azithromycin",
    #               "m_amoxicillin", "m_cefazolin","m_erythromycin","m_gentamicin","m_penicillin","m_vancomycin",
    #               "m_clindamycin","m_cefotaxime", "dur_membrane_rupture", "CRIB II Score","Total Antibiotic Days", "Antibiotic_Treatment","Antibiotic_Treatment_unfiltered", "Cohort"]

    all_columns = df.columns
    species_list = []
    # Get the list of unique species:
    for col in all_columns:
        species = col.replace("#abundance", "").replace("#PTR", "")
        if species not in species_list:
            species_list.append(species)


    df["AveragePTR"] = df.apply(lambda row: combine(row, species_list), axis=1)
    df["sample"] = df.index
    df = df.reset_index(drop=True)
    # Combine with metadata:
    metadata_df = pd.read_csv(metadata_path)

    df = metadata_df.merge(df, how="right", on="sample")

    df = df[["AveragePTR"] + relevant_coluns]

    df.to_csv(ptr_combined, index=False)


    # ### 5 - Select Major Species
    # Select taxa for which we have PTR values in at least 20 samples

    # In[5]:


    in_file = intermediate_dir+"/PTR_species_filtered.csv"
    metadata_file = intermediate_dir+"/PTR_sample_filtered.csv"

    out_file_major = out_dir+"/PTR_species_filtered_metadata_major.csv"
    out_file_all = out_dir+"/PTR_species_filtered_metadata.csv"

    cutoff = 20

    # Filter taxa that have less than "cutoff" values:

    # Get list of species:
    ptr_df = pd.read_csv(in_file)
    columns = ptr_df.columns
    species_list=[]
    for col in columns:
        if "PTR" in col:
            species = col.replace("#PTR", "")
            species_list.append(species)
    print("There are {} species in the dataset.".format(len(species_list)))


    filter_major(0, out_file_all, ptr_df, species_list, metadata_file)
    merged_df = filter_major(cutoff, out_file_major, ptr_df, species_list, metadata_file)


    # In[6]:


    merged_df


    # In[7]:


    # Re-normalize the abundance (some species for which PTR was not computed have been droped)
    all_cols = merged_df.columns
    abundance_cols = []
    for x in all_cols:
        if "#abundance" in x:
            abundance_cols.append(x)

    abundance_df = merged_df[abundance_cols]
    abundance_df.index = merged_df['sample']
    
    abundance_df_t = abundance_df.T
    for col in abundance_df_t.columns:
        tmp_sum = abundance_df_t[col].sum()
        abundance_df_t[col] = abundance_df_t[col].apply(lambda x: x/tmp_sum)
    abundance_df_t
    abundance_df_t['avg'] = abundance_df_t.mean(axis=1)
    abundance_df = abundance_df_t.T
    abundance_df


    # In[8]:


    for col in abundance_cols:
        print('"{}",'.format(col))


    # In[9]:


    merged_df.index = merged_df['sample']
    for x in all_cols:
        if "#abundance" in x:
            merged_df[x] = abundance_df[x]
    merged_df.to_csv(out_dir+"/PTR_species_filtered_metadata_major.csv", index=False)
    merged_df


    # In[ ]:




