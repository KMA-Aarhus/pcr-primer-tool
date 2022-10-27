import glob

list_of_files = (glob.glob("input/*"))
list_of_rids = [file.partition("/")[2].partition("-")[0] for file in list_of_files]

rule all:
	input: expand(["output/{rid}.csv"], rid = list_of_rids)

rule separate_ncbi_file:
	input:
		"input/{rid}-Alignment.txt"
	output:
		description = temp("output/{rid}_description.txt"),
		alignment = temp("output/{rid}_alignment.txt"),
		alignment_header = temp("output/{rid}_alignment_header.txt")
	shell: """

	# Make description file with first 3 words of description + accession no
	awk '/%/ {{print $1,$2,$3, $NF}}' {input} > {output.description}

	# Make alignment file without reference
	# Requires two dots in the alignment column. Requires alignment column to be the 3rd column

	awk 'index($3, "..")' {input} > {output.alignment}

	# Make header file for alignment
	awk '{{ if($1 == "Query") print}}' {input} > {output.alignment_header}

	"""

rule create_csv_file:
	input:
		description = "output/{rid}_description.txt",
		alignment = "output/{rid}_alignment.txt",
		alignment_header = "output/{rid}_alignment_header.txt"
	params:
		max_number_of_nucleotides = 150
	output:
		"output/{rid}.csv"
	script:
		"make_csv.py"