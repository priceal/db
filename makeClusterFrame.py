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
clusterResTsv = 'clusterRes4687_cluster.tsv'

##############################################################################
df=pd.read_csv(clusterResTsv,sep='\t',names=['cluster','member'])
clusterNames = set(df['cluster'] )

#=df['cluster'].value_counts()
#cv.hist(bins=1000)
#plt.semilogy()

clusterDict = { 'cluster':[], 'rep':[], 'species':[] }
for name in clusterNames:
    path = os.path.join(fastaDirectory,name[:4].lower()+'.fasta')
    read = SeqIO.parse( path , 'fasta')
    print(name,end=' ')
    for record in read: 
        if record.id.split('|')[0]==name:
            clusterDict['cluster'].append(name)
            clusterDict['rep'].append(record.description.split('|')[-2])
            clusterDict['species'].append(record.description.split('|')[-1])
clusterDf = pd.DataFrame(clusterDict)