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

can limit pdbids to those in a filtered list if wanted.

"""

import requests
import pandas as pd
import os

# define script variables
# inputs
summaryFile = './csv/summary.csv'  # name of summary file 
maxNumber = 1000      # maximum number to download (limit to first maxNumber ids)
filterFile = '../PDB/csv/sort713.csv'  # leave '' for no filtering, omly those
# included in this file will be downloaded

#outputs
assemblyDirectory = '../DATA/db/assemblies'  # directory to contain assembly files

###############################################################################
# load in summary dataframe
summaryDf = pd.read_csv(summaryFile)
print(summaryFile,':',len(summaryDf),'total entries')
print('limiting to first',maxNumber,'entries')

# read in filterFile if defined, and create filter set
if filterFile:
    filterDf = pd.read_csv(filterFile)
    filterSet = { p.lower() for p in filterDf['pdbid'] }

'''
the following will create the directories for structure of the preferred 
assemblies. and download all from the PDB.
'''
os.makedirs(assemblyDirectory,exist_ok=True)
print('downloading data files to', assemblyDirectory)
skippedPdbids = []
downloadCount = 0
for entry in summaryDf.itertuples():
    if downloadCount == maxNumber:
        break
    elif entry.pdbid not in filterSet:
        skippedPdbids.append(entry.pdbid)
    else:
        fileName = entry.pdbid + '-assembly' + str(entry.assembly_id) + '.cif'
        url = 'https://files.rcsb.org/download/' + fileName
        download = requests.get(url)
        with open( os.path.join(assemblyDirectory,fileName), 'w' ) as f:
            f.write( download.text )
        downloadCount += 1
        print(entry.pdbid,end=' ')
print(f'\ndownload {downloadCount} PDB assembly files')
print(f'\nskipped {len(skippedPdbids)} pdbids:',end=' ')
print(' '.join(skippedPdbids))


        
        