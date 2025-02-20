#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 08:28:07 2025

perform an analysis of header information in a directory of mmCIF files.
for current purposes, assume they are ASU files, not assemblies

@author: allen
"""
from Bio import SeqIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict    # parses data in mmCIF files
import pandas as pd
import os

# map CIF tokens onto dictionary keys. 
# note: not all tokens of interest appear in all entries.
# in definitions below, key is name of field (column) in database (dataframe), and value is the CIF token.
tokens = {'pdbid': '_entry.id',\
          'title': '_struct.title',\
          'method': '_exptl.method',\
          'date': '_pdbx_database_status.recvd_initial_deposition_date',\
          'keywords': '_struct_keywords.pdbx_keywords',\
          'text': '_struct_keywords.text',\
          'related': '_pdbx_database_related.db_id',\
          'entity type': '_entity.type',\
          'entity descr': '_entity.pdbx_description',\
          'entity MW': '_entity.formula_weight',\
          'entity number': '_entity.pdbx_number_of_molecules',\
          'entity mut': '_entity.pdbx_mutation',\
          'poly type': '_entity_poly.type',\
#          'poly seq': '_entity_poly.pdbx_seq_one_letter_code_can',\
          'poly chain': '_entity_poly.pdbx_strand_id',\
          'assembly Oligo': '_pdbx_struct_assembly.oligomeric_details',\
          'assembly Count': '_pdbx_struct_assembly.oligomeric_count'
          }
    
# define data directory
pdbDirectory = 'pdb'
#fastaDirectory = 'fasta'

dirList = os.listdir(pdbDirectory)
dirList.sort()
pdbCodes = {}
for s in dirList:
    pdbCodes[s[:4]]=s
print(len(pdbCodes),'codes')

# create the dictionaries of the values I want from the headers 
dataDict = { k:[] for k in tokens.keys() }    # holds the values for each entry
lengthDict = { k:[] for k in tokens.keys() }  # holds counts of each value (length of list)

# parse mmCIF file, and map values to dictionaries
print('processing data directory:', pdbDirectory)
for code,file in pdbCodes.items():
    print(code,end=' ')
    mmCifDict = MMCIF2Dict(os.path.join(pdbDirectory,file))
 #   sequence=SeqIO.parse(os.path.join(fastaDirectory,code+'.fasta'), "fasta")
    for k,v in tokens.items():        # iterate over asu token/value pairs
        try:
            value = mmCifDict[v]
        except:
            value = []
        dataDict[k].append(value)
        lengthDict[k].append(len(value))
#    print(sequence)
    
dataDf=pd.DataFrame(dataDict)
lengthDf=pd.DataFrame(lengthDict)








