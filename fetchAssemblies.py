#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 10:40:36 2025

@author: allen

database:   Protein Data Bank 
            https://files.rcsb.org

downloads preferred assembly structure for a list of PDB ids. 
Needs the summary dataframe saved as a csv created by 
'downloadPDBeSummary.py', as it reads the PDB ids and the preferred 
assembly label from that dataframe.

"""

import requests
import pandas as pd
import os

# define script variables
summaryFile = 'summary.csv'  # name of summary file 
assemblyDirectory = '../DATA/db/assemblies'   # directory to contain assembly files
batchSize = 100            # batch size for download (must be <1000)
maxNumber = 20      # maximum number to download (limit to first maxNumber ids)

###############################################################################
# load in dataframe
summaryDf = pd.read_csv(summaryFile)
print(summaryFile,':',len(summaryDf),'total entries')
print('limiting to first',maxNumber,'entries')
'''
the following will create the directories for structure of the preferred 
assemblies. and download all from the PDB.
'''
os.makedirs(assemblyDirectory,exist_ok=True)
print('downloading data files to', assemblyDirectory)
for i in summaryDf.index[:maxNumber]:
    code=summaryDf.at[i,'pdbid']
    assembly=summaryDf.at[i,'assembly_id']
    fileName = code + '-assembly' + str(assembly) + '.cif'
    url = 'https://files.rcsb.org/download/' + fileName
    download = requests.get(url)
    with open( os.path.join(assemblyDirectory,fileName), 'w' ) as f:
        f.write( download.text )

print('download complete.')
        
        