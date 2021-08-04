#! /bin/python

"""
This script is used to merge multiple vcf files into a single csv file.
The script is written in Python 3.7.3.
Written by Margaret Wanjiku on 07 Dec 2019.

Usage: python vcf_merger.py <in_folder> <out_folder>
    in_folder: Path to folder containing vcf files to be merged.
    out_folder: Path to folder where csv file is written.
    
NB: Do not use absolute path when specifying the paths.
"""

# Import required modules
import os
import sys
import pandas as pd
import allel

# Check of all arguments have been provided
if len(sys.argv) == 3:
    # Input and output folder paths
    in_path = os.getcwd() + "/" + sys.argv[1]
    out_path = os.getcwd() + "/" + sys.argv[2]
    
    # Check if the paths exist
    if os.path.isdir(in_path):
        pass
    else:
        print("The input folder does not exist!!!")
        print(__doc__)
        exit()
        
    if os.path.isdir(out_path):
        pass
    else:
        print("The output folder does not exist!!!")
        print(__doc__)
        exit()
    
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
    
    #Write output to a csv file
    print("Writing to csv file " + out_path + 'merged_vcf.csv')
    merged_vcf.to_csv(out_path + 'merged_vcf.csv', index = False)

    print("Merged successfully !!!")
    
else:
    print(__doc__)
