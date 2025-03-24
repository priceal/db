# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 10:43:12 2025

@author: priceal
"""

def is_dna( seq ):
    result = True
    for c in seq:
        result*=(c in {'A','G','C','T'})
    return bool(result)