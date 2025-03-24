#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 10:40:36 2025

@author: allen

database:   Protein Data Bank 
            https://files.rcsb.org

Downloads sequence files for a list of PDB ids. The sequence data
is the fasta file from PDB--containing the sequence of the crystallized
entity in the case of x-ray strutures, with sequences of all polymer entities.

N.B. This sequence will contain all residues, even those not observed in
experimental structure.

List of PDB ids can be either a .csv with a column header 'pdbid', created 
by downloadPDBeSummary.py for example. Or it can be a white space separated
list of PDB ids.

"""

import requests
import json
import pandas as pd
from itertools import batched
import os
from Bio import SeqIO

# define script variables            
fastaDirectory = '../DATA/db/fasta'            # directory to contain fasta files
maxNumber = 20      # maximum number to download (limit to first maxNumber ids)


###############################################################################
def is_dna( seq ):
    for c in seq:
        if c not in 'ACGT':
            return False
    return True

def is_protein( seq ):
    for c in seq:
        if c not in 'ARNDCEQGHILKMFPSTWYVX':
            return False
    return True
###############################################################################

# create list of files paths
dirList = os.listdir(fastaDirectory)
pathList = [ os.path.join(fastaDirectory,f) for f in dirList ]
pdbids = [ p[:4] for p in dirList ]
dataDf = pd.DataFrame({'pdbid':pdbids, 
                       'path':pathList,
                       'dna':[[]]*len(pdbids),
                       'protein':[[]]*len(pdbids),
                       'duplex': [[]]*len(pdbids)}
                      )
   
'''
the following will create the directory for sequence files and download the 
fasta files. Here, all are downloaded from the PDB.
'''

for i in dataDf.index:
    record = SeqIO.parse( dataDf.at[i,'path'], 'fasta' )
    print( dataDf.at[i,'path'] )
    d =[]; p = []; entryDNAs=[]
    for entry  in record:
        print('\t',entry.id, is_dna( entry.seq),len(entry.seq),end=' ')
        if is_dna( entry.seq ):
            d.append(len(entry.seq))
            entryDNAs.append(entry)
            print( entry.seq, entry.reverse_complement().seq)
        elif is_protein( entry.seq ):
            p.append(len(entry.seq))
            print('')
        else:
            print('unrecognized residue')
        
    dataDf.at[i,'dna']=d
    dataDf.at[i,'protein']=p
    if len( entryDNAs ) == 1:
        duplex= (entryDNAs[0].seq==entryDNAs[0].reverse_complement().seq)
    elif len( entryDNAs ) ==2:
        duplex= (entryDNAs[0].seq==entryDNAs[1].reverse_complement().seq)
    dataDf.at[i,'duplex'] = duplex
    
            
    # the fasta file - n.b. need upper case for these files!

        
        