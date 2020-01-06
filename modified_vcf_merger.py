# Import required modules
import os
import sys
import pandas as pd
import allel
import numpy as np

# Input and output folder paths
in_path = os.getcwd() + "/" 
out_path = os.getcwd() + "/"

#List files in that directory
path = in_path # path to vcf files to be merged
files = os.listdir(path) # list files in a directory
#print(files)

#Select only vcf files 
vcf_files = []
for file in files:
    if file.split(".")[-1].casefold() == "vcf":
        vcf_files.append(file)
        #print(file)
        
vcf_files.sort()
#print(vcf_files)


### START THE MERGING PROCESS
print("Started merging...")

# Read the first two files
vcf_1 = allel.vcf_to_dataframe(path + "/" + str(vcf_files[0]), fields=['CHROM', 'POS', 'REF', 'ALT'], alt_number=1)
vcf_2 = allel.vcf_to_dataframe(path + "/" + str(vcf_files[1]), fields=['CHROM', 'POS', 'REF', 'ALT'], alt_number=1)

#Dropping rows with insertions e.g ACCT, AT
vcf_1 = vcf_1[vcf_1['ALT'].str.len().lt(2)]
vcf_2 = vcf_2[vcf_2['ALT'].str.len().lt(2)]

print("Merging " + str(vcf_files[0]))
print("Merging " + str(vcf_files[1]))

# Merge the first two vcf files
merged_vcf = pd.merge(vcf_1, vcf_2, how = 'outer', left_on = ['CHROM', 'POS', 'REF'], right_on = ['CHROM', 'POS', 'REF'])

#Rename the columns using a dictionary
merged_vcf.rename(columns={
        'ALT_x': str(vcf_files[0].split('_')[-1].split('.')[0]), 'ALT_y': str(vcf_files[1].split('_')[-1].split('.')[0])
    }, inplace=True)

#Loop through the other vcf files
for vcf in vcf_files[2:]:
    # Read the vcf file
    tmp_vcf = allel.vcf_to_dataframe(path + "/" + str(vcf), fields=['CHROM', 'POS', 'REF', 'ALT'], alt_number=1)

    #Drop insertions
    tmp_vcf = tmp_vcf[tmp_vcf['ALT'].str.len().lt(2)]
    print("Merging " + str(vcf))

    # Merge the vcf file to the master file
    merged_vcf = pd.merge(merged_vcf, tmp_vcf, how = 'outer', left_on = ['CHROM', 'POS', 'REF'], right_on = ['CHROM', 'POS', 'REF'])

    # Rename the added column
    merged_vcf.rename(columns={
            'ALT': str(vcf.split('_')[-1].split('.')[0])
        }, inplace=True)
    
#Sort the merged vcf file
merged_vcf = merged_vcf.sort_values(['CHROM','POS'], ascending=[True,True])

#Replace all the NaN values with the string, NAN
merged_vcf.replace(np.nan, 'NAN', inplace = True)

#Replace all the NAN strings in each isolate with the string in the corresponding reference column
lst = []
for col in merged_vcf.columns:
    lst.append(col)

for i in lst[3:]:
    for j in merged_vcf.index:
        if merged_vcf.at[j, i] == "NAN":
            merged_vcf.at[j, i] = merged_vcf.at[j, 'REF']
            
        else:
            pass


merged_vcf.T
transposed_merged_vcf_file = merged_vcf.T
#new_file.csv = transposed_merged_vcf_file.to_csv

#Write output to a csv file
print("Writing to csv file " + out_path + 'new_file.csv')
transposed_merged_vcf_file.to_csv(out_path + 'new_file.csv')


print("Merged successfully !!!")
