AZ scripts
==============

The script cached_entrez.py contains funcion `efetch_multiple` that retrieves
fasta and gb info from the nucleotide Genbank database.
It accepts a query, saves matched sequences to files and returns the dirpath.
If you call it the second time, it will just return the dirpath without fetching
the sequencing once more.
See fetch_proteomes.py for the usage.

The folder othomcl_software contains all OrthoMCL stuff, including binaries which should work for all
UNIX system.

The orthomcl.config is properly corrected configuration file for OrthoMCL.