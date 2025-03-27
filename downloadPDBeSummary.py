#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 10:40:36 2025

@author: allen

database:   Protein Data Bank in Europe
            https://www.ebi.ac.uk/pdbe/
            
downloads summary information for a list of PDB ids. the summary information 
is downloaded from the PDBe and is used to create a summary dataframe which 
can be saved in a csv file. 

"""

import requests
import json
import pandas as pd
from itertools import batched
import os

# define script variables
# inputs
pdbCodeFile = 'csv/pdbListAll.txt'              # file containing pdb ids
batchSize = 100            # batch size for download (must be <1000)
maxNumber = 1000      # maximum number to download (limit to first maxNumber ids)

# outputs
summaryFile = 'csv/summary.csv'  # name of summary file to save, set to '' if not wanted

###############################################################################
# load in pdb ids
# use below if csv file with one column labeled 'pdbid'
if os.path.splitext(pdbCodeFile)[-1] == '.csv':  
    df = pd.read_csv(pdbCodeFile)
    pdbCodes=list(df['pdbid'])
else:    # a simple whitespace separated list of ids
    with open(pdbCodeFile) as f:
        fileRead=f.read()
    pdbCodes = fileRead.strip().split()

# limit to maxNumber
print('read in',len(pdbCodes),'pdb ids\n')
pdbCodes = pdbCodes[:maxNumber]
print('limiting to first',len(pdbCodes),'pdb ids\n')

# download multiple entries from PDBe, taking care that PDBe only allows batches 
# up to 1000 at a time
urlPrefix = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/'
pdbCodesBatched = batched( pdbCodes, batchSize ) # creates list of batches 
reportDict = {}   # to store downloads
pdbCodes = [c.lower() for c in pdbCodes] # need lower case for PDBe
print('downloading summary data from PDBe')
for batch in pdbCodesBatched:
    print('downloading batch...', end='')
    # data=list of comma separated ids
    report=requests.post(urlPrefix,data=','.join(batch)) 
    # convert json format to dictionary and append to dictionary
    reportDict.update( json.loads(report.text) ) 
    print(len(reportDict),'total entries')
print(len(reportDict),'entries downloaded')

'''
the following creates the summary dataframe from downloaded dictionary. First
are lists of various keys used in the dictionaries.

# keys in downloaded dictionary 'reportDict' (not used but listed for completeness)
entryKeys = ['title', 'processing_site', 'deposition_site', 'deposition_date', \
             'release_date', 'revision_date', 'experimental_method_class', \
             'experimental_method', 'split_entry', 'related_structures', \
             'entry_authors', 'number_of_entities', 'assemblies']
'''
# keys with numerical, string or list values
simpleKeys = ['title', 'deposition_date', 'experimental_method']

# keys of the sub-dictionaries that are the values associated with keys 
# 'number_of_entities' and 'assemblies'
entityKeys = ['water', 'polypeptide', 'dna', 'rna', 'dna/rna', 'sugar', \
              'ligand', 'carbohydrate_polymer', 'other']
assemblyKeys = ['assembly_id', 'name', 'form']

# keys (columns) of summary dictionary (dataframe). 'assemblies' will hold the
# number of assemblies for the entry
dataKeys = ['pdbid'] + simpleKeys + entityKeys + ['assemblies'] + assemblyKeys

# creation of the summary dictionary, 'dataDict'. Must use entry[0] to extract
# dictionary as it is surrounded by [] when returned by json.loads()
dataDict = { k:[] for k in dataKeys }   # blank dictionary of lists
print('creating summary dataframe...')
for pdbid,entry in reportDict.items():
    dataDict['pdbid'].append(pdbid)
    for k in simpleKeys:
        dataDict[k].append(entry[0][k])
    for k in entityKeys:
        dataDict[k].append(entry[0]['number_of_entities'][k])
    dataDict['assemblies'].append(len(entry[0]['assemblies']))

    # now go through the assemblies and extract data from preferred assembly: 
    # the PDBe defines the preferred assembly as the smallest assembly 
    # containing all polymeric entities.
    for d in entry[0]['assemblies']:
        if d['preferred']:    # this is True only for preferred assembly
            for ak in assemblyKeys:
                dataDict[ak].append(d[ak])
dataDf = pd.DataFrame(dataDict)  # create the dataframe

# print for review and save
print(dataDf)
print( dataDf.describe() )
if summaryFile: dataDf.to_csv(summaryFile) # save if requested


        