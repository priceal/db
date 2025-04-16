#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 15:30:48 2025

@author: allen
"""

import os
import pandas as pd
import pylab as plt
from Bio import SeqIO

# inputs
fastaDirectory = '../DATA/db/fasta'


df=pd.read_csv('clusterRes4687_cluster.tsv',sep='\t',names=['cluster','member'])
cv=df['cluster'].value_counts()
cv.hist(bins=1000)
plt.semilogy()

clusterNames = set(df['cluster'] )
clusterDict = { 'cluster':[], 'rep':[], 'species':[] }
for name in clusterNames:
    path = os.path.join(fastaDirectory,name[:4].lower()+'.fasta')
    print(f'\n{name = }')
    read = SeqIO.parse( path , 'fasta')
    for record in read: 
        if record.id.split('|')[0]==name:
            print(record.description)
            clusterDict['cluster'].append(name)
            clusterDict['rep'].append(record.description.split('|')[-2])
            clusterDict['species'].append(record.description.split('|')[-1])
clusterDf = pd.DataFrame(clusterDict)