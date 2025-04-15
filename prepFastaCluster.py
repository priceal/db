#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 10:58:44 2025

@author: allen
"""

import os
import random
from Bio import SeqIO
from Bio.PDB import is_aa

fastaDirectory = '../DATA/db/fasta'
sampleSize = 100
logFile = 'prepFastaCluster_test_20250414.log'


onlyAminoAcids = "RNDEQHILKMFPSWYV"

###############################################################################
def is_protein( chain ):
    '''
    determines if at least 5 residues in chain are canonical AA residue types 

    Args:
        chain (TYPE): DESCRIPTION.

    Returns:
        bool: DESCRIPTION.

    '''
    count = 0
    for r in chain:
        if r in onlyAminoAcids: count += 1
            
    return count > 1

###############################################################################

'''
###############################################################################
########################    main   ##########################################
###############################################################################
'''

listDir = os.listdir(fastaDirectory)
sampleList = random.sample(listDir,sampleSize)

with open(logFile,'w') as f:

    # load in all fasta files in list, and create a single fasta file
    # with all sequences and save.
    outputList= []
    for file in sampleList:
        read = SeqIO.parse(os.path.join(fastaDirectory,file), 'fasta' )
        f.write('-----file: '+file+'\n')
        for record in read:
            f.write(f'{record.id}\nlength = {len(record)}\n')    
            f.write(repr(record.seq))
            if is_protein( str(record.seq) ):
                outputList.append(record)
                f.write(f'    ....protein! count at {len(outputList)}\n')

# now run mmseqs2
SeqIO.write(outputList, "temp.fasta", "fasta")
