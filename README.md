AZ scripts
==============

The `orthomcl.config` is a properly corrected configuration file for OrthoMCL.

The script `cached_entrez.py` contains the function `efetch_multiple` that retrieves
fasta and gb info from the nucleotide Genbank database.
It accepts a species name, then saves matched sequences to files into the `cache`
directory and returns the paths.
If you call it the second time, it will just return the paths without fetching
the sequencing once more.

The `fetch_proteomes.py` script utilizes `efetch_multiple` to retrieve genomes,
then iterates CDS annotations and writes OrthoMCL-compliant translations into
the `../data/<species name>_proteomes` directory.

The `blastp_proteomes.py` script performs filters fasta with orthomclFilterFasta,
and performs blastp all versus all, produsing blasted.xml file for OrthoMCL.


Requirements:
pip install mysql-connector-python