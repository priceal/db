#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 10:40:36 2025

@author: allen

database:   Protein Data Bank 
            https://files.rcsb.org

Downloads sequence files for a list of PDB ids. The sequence data
is the fasta file from PDB--containing the sequence of the crystallized
entity in the case of x-ray strutures, with sequences of all polymer entities.

N.B. This sequence will contain all residues, even those not observed in
experimental structure.

List of PDB ids can be either a .csv with a column header 'pdbid', created 
by downloadPDBeSummary.py for example. Or it can be a white space separated
list of PDB ids.

"""

import requests
import pandas as pd
import os

# inputs
pdbCodeFile = './csv/filtered_20250402.csv'              # file containing pdb ids

#outputs
fastaDirectory = '../DATA/db/fasta'            # directory to contain fasta files
maxNumber = 10000      # maximum number to download (limit to first maxNumber ids)

###############################################################################
# load in pdb ids
# use below if csv file with one column labeled 'pdbid'
if os.path.splitext(pdbCodeFile)[-1] == '.csv':  
    df = pd.read_csv(pdbCodeFile)
    pdbCodes=list(df['pdbid'])
else:    # or below if a simple whitespace separated list of ids
    with open(pdbCodeFile) as f:
        fileRead=f.read()
    pdbCodes = fileRead.strip().split()
    
# limit to maxNumber
print('read in',len(pdbCodes),'pdb ids\n')
pdbCodes = pdbCodes[:maxNumber]
print('limiting to first',len(pdbCodes),'pdb ids\n')

'''
the following will create the directory for sequence files and download the 
fasta files. Here, all are downloaded from the PDB.
'''
os.makedirs(fastaDirectory,exist_ok=True)
print('downloading fasta files to', fastaDirectory)
# n.b. convention is to ensure lower case for all file names
# however, PDB  uses upper case for fasta API
for code in pdbCodes:
    savePath = os.path.join(fastaDirectory,code.lower()+'.fasta')
    if os.path.exists(savePath):
        print('*'+code+'*',end=' ') # indicate skipped
    else:
        url = 'https://www.rcsb.org/fasta/entry/' + code.upper()
        try:
            download = requests.get(url)
            with open( savePath, 'w' ) as f:
                f.write( download.text )
            print(code,end=' ')
        except Exception as e:
            print(f"Error fetching {code}: {e}")
            
print('\ndownloads completed.')
        
        
        