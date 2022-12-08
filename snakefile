# snakemake --profile configs/profile

query_sequence = config["query_sequence"]
query_organism = config["query_organism"]
file_id = query_organism.replace(" ", "_")


rule all:
	input: ["output/" + file_id + "_all_alignments.csv"]


rule download_database:
	output: "database/" + file_id + ".fasta"
	params:
		query = query_organism
	conda: "configs/envs/ncbi.yaml"
	script: "scripts/database_download.py"

rule blast:
	input: 
		database = "database/" + file_id + ".fasta"
	output:
		alignment_file = temp(file_id + "_alignment.txt"),
		description_file = temp(file_id + "_description.tsv"),
		query_seq = temp(file_id + "query.txt")
	params:
		database_title = file_id
	conda: "configs/envs/ncbi.yaml"
	shell: """

	# Prepare database for BLAST
	makeblastdb -in {input.database} -dbtype nucl -parse_seqids -out {params.database_title}

	# Use blastn -help for list of BLAST output formats and format specifier options

	echo {query_sequence} > {output.query_seq}

	# BLAST alignment file
	blastn -query {output.query_seq} -db {params.database_title} -outfmt 3 -line_length 1000 -max_target_seqs 10000 -out {output.alignment_file} 

	# BLAST description. Creates tsv with 3 columns: Subject accession.version, Subject Title & Percentage of identical matches
	blastn -query {output.query_seq} -db {params.database_title} -outfmt "6 saccver stitle pident" -max_target_seqs 1000 -out {output.description_file}

	"""


rule separate_ncbi_file:
	input:
		file_id + "_alignment.txt"
	output:
		alignment = temp(file_id + "_alignments.txt"),
		alignment_header = temp(file_id + "_alignment_header.txt")
	shell: """

	# Make alignment file without reference
	# Requires alignment column to be the 3rd column. Requires four dots in the alignment column. 
	awk 'index($3, "....")' {input} > {output.alignment}

	# Make header file for alignment
	awk 'index($1, "Query_")' {input} > {output.alignment_header}

	"""

rule create_csv_files:
	input:
		description = file_id + "_description.tsv",
		alignment = file_id + "_alignments.txt",
		alignment_header = file_id + "_alignment_header.txt"
	output:
		"output/" + file_id + "_all_alignments.csv",
		"output/" + file_id + "_alignments_with_nucleotide_changes.csv",
		"output/" + file_id + "_nucleotide_changes_table.csv"
	script:
		"scripts/make_csv.py"

rule clean:
	shell: 'mv {file_id}* /database'