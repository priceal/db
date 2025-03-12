
Tools for fetching and parsing data from bioinformatics databases -
mainly structural

###############################################################################
downloadPDBeSummary.py
###############################################################################

database:   Protein Data Bank in Europe
            https://www.ebi.ac.uk/pdbe/
            
Downloads summary information for a list of PDB ids. the summary information 
is downloaded from the PDBe and is used to create a summary dataframe which 
can be saved in a csv file. 

###############################################################################
fetchAssemblies.py
###############################################################################

database:   Protein Data Bank 
            https://files.rcsb.org

Downloads preferred assembly structure for a list of PDB ids. 
Needs the summary dataframe saved as a csv created by 
'downloadPDBeSummary.py', as it reads the PDB ids and the preferred 
assembly label from that dataframe.

###############################################################################
fetchFasta.py
###############################################################################

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
