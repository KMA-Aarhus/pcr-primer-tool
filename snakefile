import glob

list_of_files = (glob.glob("input/*"))
list_of_rids = [file.partition("/")[2].partition("-")[0] for file in list_of_files]

rule all:
	input: expand(["output/{rid}_all_alignments.csv", "database.fasta", "alignment.txt"], rid = list_of_rids)


rule download_database:
	output: "database.fasta"
	params:
		query = "Monkeypox"
	script:
		"database_download.py"

rule blast:
	input: 
		database = "database.fasta"
		query = "query.txt"
	params:
		database_title = "database"
	output:
		alignment_file = "alignment.txt"
		description_file = "description.txt"
	shell: """

	# Prepare database for BLAST
	makeblastdb -in {input.database} -dbtype nucl -parse_seqids -out {params.database_title}

	# BLAST create alignment file
	blastn -query {input.query} -db {params.database_title} -outfmt 3 -line_length 1000 -max_target_seqs 10000 -out {output.alignment_file} 

	#BLAST description
	blastn -query {input.query} -db {params.database_title} -outfmt "6 sacc stitle pident" -max_target_seqs 1000 -out {output.description_file}

	"""


rule separate_ncbi_file:
	input:
		"input/{rid}-Alignment.txt"
	output:
		description = temp("output/{rid}_description.txt"),
		alignment = temp("output/{rid}_alignment.txt"),
		alignment_header = temp("output/{rid}_alignment_header.txt")
	shell: """

	# Make description file with full description, percent identity and accession no
	awk '/%/ {{cols="\t"$(NF-2)"\t"$NF; NF-=11; print $0 cols}}' {input} > {output.description}

	# Make alignment file without reference
	# Requires four dots in the alignment column. Requires alignment column to be the 3rd column
	awk 'index($3, "....")' {input} > {output.alignment}

	# Make header file for alignment
	awk '{{ if($1 == "Query") print}}' {input} > {output.alignment_header}

	"""

rule create_csv_files:
	input:
		description = "output/{rid}_description.txt",
		alignment = "output/{rid}_alignment.txt",
		alignment_header = "output/{rid}_alignment_header.txt"
	params:
		max_number_of_nucleotides = 150
	output:
		"output/{rid}_all_alignments.csv",
		"output/{rid}_alignments_with_nucleotide_changes.csv",
		"output/{rid}_nucleotide_changes_table.csv"
	script:
		"make_csv.py"