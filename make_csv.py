import pandas as pd
import numpy as np

#input: txt file downloaded from ncbi, after alignment view has been changed to dots for identity

#when making the workflow, remember to remove intermediary files

#########################
# CREATE DESCRIPTION DF #
#########################

# Read file
desc = pd.read_csv(snakemake.input[0], header = None, sep = " ")

# Combine columns
desc[4] = desc[[0,1,2]].agg(' '.join, axis=1)
desc = desc.drop(desc.iloc[:, 0:3],axis = 1)

# Set header names
desc.columns = ["Query", "Description"]

assert desc["Query"].is_unique, "Description dataframe contains duplicate queries. Please remove these before continuing."
#could maybe make some code that creates dataframe with unique values. Not sure if this will ever be a problem.

#################
# CREATE HEADER #
#################

def concatenate_reference_strings(header):
	reference = ""
	for row in range(len(header)):
		reference += header.iloc[row,2]
	return reference

def create_dataframe_from_reference_string(reference):
	df = pd.DataFrame()
	df[0] = [reference]
	df = df[0].str.split(pat = "", n = len(reference), expand = True)
	return df

def concatenate_header_columns(header, reference_df):
	reference_df.insert(len(reference_df.columns), len(reference_df.columns), len(reference_df.columns)-1)
	header = pd.concat([header.iloc[0:1, 0:2], reference_df.iloc[:, 1:]], axis = 1)
	return header

def create_header(headerfile):
	header = pd.read_csv(headerfile, header = None)
	header = header[0].str.split(expand = True)
	reference = concatenate_reference_strings(header)
	reference_df = create_dataframe_from_reference_string(reference)
	header = concatenate_header_columns(header, reference_df)
	return header.iloc[0]

header = create_header(snakemake.input[2])

#######################
# CREATE ALIGNMENT DF #
#######################

def split_query_column_from_alignment(alignment):
	return alignment[0].str.split(n = 1, pat = "  ", expand = True)

def create_array_of_groups(query_alignment):
	total_rows = len(query_alignment[0])
	unique_rows = len(query_alignment[0].unique())
	n_groups = int(total_rows/unique_rows)
	list_of_groups = np.arange(1, n_groups + 1)
	groups = np.repeat(list_of_groups, total_rows)
	return groups

def pivot_rows_to_columns(query_alignment, group_by_column):
	pivoted_alignment = query_alignment.pivot(index = 0, columns = group_by_column, values = 1)
	pivoted_alignment.reset_index(inplace = True)
	return pivoted_alignment

def create_alignment_df(alignment_file, header):
	alignment = pd.read_csv(alignment_file, header = None)
	# THIS FUNCTION CURRENTLY ONLY WORKS WHEN THERE ARE NO DUPLICATES IN EITHER DESC OR ALIGN
	assert len(desc) == len(alignment[0].unique()), f'Expected {len(desc)} alignments, got {len(alignment[0].unique())}'

	query_alignment = split_query_column_from_alignment(alignment)
	groups = create_array_of_groups(query_alignment)
	query_alignment['groups'] = groups.tolist()
	pivoted_alignment = pivot_rows_to_columns(query_alignment, 'groups')

	final_alignment = pivoted_alignment[0]

	reference_length = header.iloc[-1]
	n_nucleotides = snakemake.params[0]
	n_groups = groups.max()

	for i in range(1, n_groups+1):
		df1 = pivoted_alignment[i].str.split(n = 1, pat = "  ", expand = True)

		start_pos = df1[0]
		if i == 1: # if this is the first column: add start position to final dataframe
			final_alignment = pd.concat([final_alignment, start_pos], axis = 1)

		if i < n_groups	or reference_length % n_nucleotides == 0: # if not last column or last column contains max number of nucleotides
			n = n_nucleotides + 1
		else: # number of nucleotides in last column varies
			n = reference_length % n_nucleotides + 1

		df2 = df1[1].str.split(pat = "", n = n, expand = True)

		if i < n_groups:# if not last column: don't include end position
			final_alignment = pd.concat([final_alignment, df2.iloc[:, 1:-1]], axis = 1)
		else: # include end position in final dataframe
			final_alignment = pd.concat([final_alignment, df2.iloc[:, 1:]], axis = 1)

	final_alignment.columns = header

	return final_alignment

alignment = create_alignment_df(snakemake.input[1], header)

####################
# MERGE DATAFRAMES #
####################

df = desc.merge(alignment, how = 'inner')

######################
# CREATE OUTPUT FILE #
######################

df.to_csv(snakemake.output[0])


