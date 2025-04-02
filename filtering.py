#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  2 13:01:42 2025

@author: allen
"""
import pandas as pd

dataDf = pd.read_csv('csv/fastaSummaries_20250402.csv',
            converters = {
                'dna': lambda x: eval(x),
                'U': lambda x: eval(x),
                'protein': lambda x: eval(x)
                }
                     )

# add column for max protein lenght and screen 
# > 19 and < 500

dataDf['max p']=dataDf['protein'].map( lambda p: max(p) )
dataP = dataDf[ dataDf['max p'] > 19 ]
dataP2 = dataP[ dataP['max p'] < 500 ]

# screen out any w/ U containing DNAs

dataP2['max U']=dataP2['U'].map( lambda p: max(p) if p else 0 )
dataU = dataP2[ dataP2['max U'] == 0]


# screen out any structures with minimum DNA length <4
# and that have non-standard bases unrecognizes sequence symbols
dataU['min dna'] = dataU['dna'].map( lambda p : min(p) if p else 0 )


dataDNA = dataU[ dataU[ 'min dna'] > 3 ]








