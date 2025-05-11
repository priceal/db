#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  1 08:36:52 2025

@author: allen
"""
import subprocess

query_fasta = '6jma.fasta'
output_hhr = '4o3n.fasta'

cmd = [
    "hhalign",
    "-i", query_fasta,          # Input query FASTA file
    "-t", output_hhr           # Output HHR file
]

print("Running HHblits with command:", " ".join(cmd))


# Run hhblits
p=subprocess.run(cmd, check=True,capture_output=True, encoding='utf-8')
