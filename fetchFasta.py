#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 08:28:07 2025

perform an analysis of header information in a directory of mmCIF files.
for current purposes, assume they are ASU files, not assemblies

@author: allen
"""
import os
import pandas as pd
import sys
sys.path.insert(0, '/home/allen/projects/PDB')
import pdbTools as pt


    
# define data directory
pdbCsv='df.csv'
dataDirectory      = 'fasta'

if not os.path.exists( dataDirectory ):
    os.mkdir( dataDirectory )

df=pd.read_csv(pdbCsv)
pdbCodes = [ s.split('\'')[1].upper() for s in df['pdbid'] ]

for code in pdbCodes:
    print( code, end=' ')
    text = pt.fetchSequence(code)
    with open( os.path.join(dataDirectory,code.lower()+'.fasta'), 'w') as f:
        f.write(text)
















