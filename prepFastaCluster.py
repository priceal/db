#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 10:58:44 2025

loads all fasta files in a directory and sifts out the protein sequences.
these are saved in a single file that can be read by other processing 
tools such as mmseqs

note: protein testing is imperfect, as is based on sequence letters.
If the chain has at least one letter from the list of symbols that only
code for amino acids (i.e., at least on that is not ACTGUN), then it's
classified protein

better way --- analyze summary info from PDBe to determine exactly which 
chains are protein.


@author: allen
"""

import os
import random
from Bio import SeqIO

# inputs
fastaDirectory = '../DATA/db/fasta'
sampleSize = 0   # non-zero to choose a subsample of files

# outputs
outputFile = 'dbFasta1000.fasta'
logFile = 'prepFastaCluster_dbFasta1000_20250416.log'

'''
###############################################################################
########################    functions   ##########################################
###############################################################################
'''
def is_protein( chain ):
    '''
    determines if at least 1 residue in chain is canonical AA residue type 
    which does not code for a DNA/RNA base. Uses a set 'trick'

    Args:
        chain (str): DESCRIPTION.

    Returns:
        bool: DESCRIPTION.

    '''
    return bool( set(chain).intersection( set('RDEQHILKMFPSWYV') ) )

'''
###############################################################################
########################    main   ##########################################
###############################################################################
'''

listDir = os.listdir(fastaDirectory)
if sampleSize==0:
    sampleList = listDir
else:
    sampleList = random.sample(listDir,sampleSize)

with open(logFile,'w') as f:

    # load in all fasta files in list, and create a single fasta file
    # with all sequences and save.
    outputList= []; count=0
    for file in sampleList:
        count+=1; print(count,end=' ')
        read = SeqIO.parse(os.path.join(fastaDirectory,file), 'fasta' )
        f.write('\n-----file: '+file+'\n')
        for record in read:
            f.write(f'\n{record.id}\nlength = {len(record)}\n')    
            f.write(repr(record.seq))
            if is_protein( str(record.seq) ):
                # terminate id at | -- simplifies id in clustering
                record.id=record.id[:record.id.find('|') ]
                outputList.append(record)
                f.write(f'    ....protein! count at {len(outputList)}\n')
                

# now run mmseqs2
SeqIO.write(outputList, outputFile, "fasta")
