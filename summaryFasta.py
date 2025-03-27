#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 10:40:36 2025

@author: allen

analyzes the fasta files in a given directory. creates a dataframe with pdbid,
path and columns for lengths of dna and proteins. Can be used for length 
filtering.

COLUMN          DESCRIPTION
dna             list of dna chain lengths
protein         list of protein chain lengths
U               list of chain lengths of chains containing 'U'
duplex          if single dna chain is palindrome, or for two chains,
                if one is reverse complement of the other

"""

import os
import pandas as pd
from Bio import SeqIO

'''
###############################################################################
#################### FUNCTIONS ################################################
###############################################################################
'''
def is_dna( seq ):
    '''
    

    Args:
        seq (TYPE): DESCRIPTION.

    Returns:
        bool: DESCRIPTION.

    '''
    for c in seq:
        if c not in 'ACGT':
            return False
    return True

def contains_U( seq ):
    return 'U' in seq

###############################################################################
def is_protein( seq ):
    '''
    

    Args:
        seq (TYPE): DESCRIPTION.

    Returns:
        bool: DESCRIPTION.

    '''
    for c in seq:
        if c not in 'ARNDCEQGHILKMFPSTWYVX':
            return False
    return True

'''
###############################################################################
#################### MAIN #####################################################
###############################################################################
'''

# inputs    
fastaDirectory = '../DATA/db/fasta'      # directory containing fasta files
maxNumber = 1000   # maximum number to process (limit to first maxNumber files)

# outputs
outputFile = 'csv/fastaSummaries.csv'  # leave '' to not save

###############################################################################

# create dataframe containing pdbids, paths and empty columns 
dirList = os.listdir(fastaDirectory)
pathList = [ os.path.join(fastaDirectory,f) for f in dirList ]
pdbids = [ p[:4] for p in dirList ]
dataDf = pd.DataFrame({'pdbid':pdbids, 
                       'path':pathList,
                       'dna':[[]]*len(pdbids),
                       'U': [[]]*len(pdbids),
                       'protein':[[]]*len(pdbids),
                       'duplex': [[]]*len(pdbids)}
                      )
    
# loop through data frame and analyze entries
for i in dataDf.index:
    record = SeqIO.parse( dataDf.at[i,'path'], 'fasta' )
    print( dataDf.at[i,'path'] )
    d =[]; u=[]; p = []; entryDNAs=[]
    for entry  in record:
        print('\t',entry.id, is_dna( entry.seq),len(entry.seq),end=' ')
        if is_dna( entry.seq ):
            d.append(len(entry.seq))
            entryDNAs.append(entry)
            print( entry.seq, entry.reverse_complement().seq)
        elif is_protein( entry.seq ):
            p.append(len(entry.seq))
            print('')
        elif contains_U( entry.seq ):
            u.append( len(entry.seq))
        else:
            print('unrecognized residue')
    # add results to dataframe       
    dataDf.at[i,'dna']=d
    dataDf.at[i,'protein']=p
    dataDf.at[i,'U']=u
    # check for duplex dna by looking for complementary ssDNAs
    if len( entryDNAs ) == 1:
        duplex= (entryDNAs[0].seq==entryDNAs[0].reverse_complement().seq)
    elif len( entryDNAs ) ==2:
        duplex= (entryDNAs[0].seq==entryDNAs[1].reverse_complement().seq)
    dataDf.at[i,'duplex'] = duplex
    
if outputFile:
    dataDf.to_csv( outputFile, index=False )
    
    