# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 10:43:12 2025

@author: priceal
"""
'''
'filtered.csv'

the 1000 files in /data/db/fasta are screened

1. add column for max protein lenght and screen 40 , max < 500

dataDf['max p']=dataDf['protein'].map( lambda p: max(p) )
dataP = dataDf[ dataDf['max p'] >40 ]
dataP2 = dataP[ dataP['max p'] <500 ]

2. add column for number of dna entities, and screen for number > 0

dataP2['number d']=dataP2['dna'].map( lambda d: len(d) )
dataD = dataP2[ dataP2['number d'] >0 ]

3. add column for min dna length, and screen for 3 < min < 50
dataD['min d'] = dataD['dna'].map( lambda d: min(d) )
dataD2 = dataD[ dataD['min d']<50 ]
dataD3 = dataD2[ dataD2['min d']>3 ]

4. require number of dnas < 3
dataD4 = dataD3[ dataD3['number d']<3 ]

4. save dataframe
dataD4.to_csv('filtered.csv')

'''