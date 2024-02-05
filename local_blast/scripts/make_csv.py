import pandas as pd
import numpy as np
import os

##################
# DESCRIPTION DF #
##################

def create_description(desc_file):
	"""
	Creates description df from a tsv file with query id, description and percent nucleotide identity columns.
	Requires query id's to be unique.
	"""
	desc = pd.read_csv(desc_file, header = None, sep = "\t")

	desc.columns = ["Query_1", "Description", "Percent_Identity"]

	desc = desc.drop_duplicates(subset =['Query_1'])
	
	return desc


##########
# HEADER #
##########

def create_header(header_file):
	"""Takes a header txt file. Creates column names for alignment dataframe"""
	header = pd.read_csv(header_file, header = None, delim_whitespace=True)
	
	# Make nucleotides uppercase
	align_col = header[2].str.upper()

	# Create one column per nucleotide
	align_cols = align_col.str.split(pat = "", expand = True)

	# Concatenate 'Query', start position, nucleotides and end position columns
	header = pd.concat([header.iloc[:,0:2], align_cols.iloc[:, 1:-1], header.iloc[:,3:]], axis = 1)

	# Get length of query sequence
	query_length = int(align_col.str.len())

	return header.iloc[0], query_length

################
# ALIGNMENT DF #
################

def create_alignment_df(alignment_file, header_file):
	alignment = pd.read_csv(alignment_file, header = None)

	header, query_length = create_header(header_file)

	# Extract Query id column
	ids = alignment[0].str.split(n = 1, expand = True)

	# Extract start position column
	start_pos = ids[1].str.split(n = 1, pat = "  ", expand = True)

	# Remove extra whitespaces due to differences in start position str length
	col = start_pos[0].str.len()
	start_pos['n_whitespaces'] = col.max() - col
	alignments = start_pos.apply(lambda x: x[1][x['n_whitespaces']:], 1)

	# Split remaining columns. One column per nucleotide position.
	alignments = alignments.str.split(pat = "", n = query_length + 1, expand = True)

	# Concatenate Query id, start position, nucleotides and end position columns
	alignment = pd.concat([ids[0], start_pos[0], alignments.iloc[:, 1:]], axis = 1)

	# Add column names
	alignment.columns = header

	return alignment

##########################################
# ERROR HANDLING IN CASE OF EMPTY SEARCH #
##########################################

def write_error_message_to_file(output_file, error_message):
	with open(output_file, 'w') as f:
		f.write(error_message)

if os.path.getsize(snakemake.input[0]) == 0:
	message = "Input files are empty. This means no hits were found in the BLAST search. Try again with a different database."

	for i in range(3):
		write_error_message_to_file(snakemake.output[i], message)

##########################
# CREATE FINAL DATAFRAME #
##########################

else:
	desc = create_description(snakemake.input[0])
	alignment = create_alignment_df(snakemake.input[1], snakemake.input[2])
	df = alignment.merge(desc, how = 'left')

#######################
# CREATE OUTPUT FILES #
#######################

# Table containing all alignments
	df.to_csv(snakemake.output[0])

# Table containing aligments with nucleotide changes
	df_nucleotide_changes = df[df["Percent_Identity"] != 100]
	df_nucleotide_changes.to_csv(snakemake.output[1])

# Table with counts of nucleotide changes
	df_count_nucleotide_changes = df.iloc[:,2:-3].apply(pd.value_counts).fillna(0).filter(regex = '[a-zA-Z]', axis=0).convert_dtypes()
	df_count_nucleotide_changes.to_csv(snakemake.output[2])