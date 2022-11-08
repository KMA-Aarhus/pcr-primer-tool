import pandas as pd
import numpy as np

#########################
# DESCRIPTION FUNCTIONS #
#########################

def create_description(desc_file):
	desc = pd.read_csv(desc_file, header = None, sep = " ")

	# Combine description columns
	desc[0] = desc[[0,1,2]].agg(' '.join, axis=1)
	desc = desc.drop(desc.iloc[:, 1:3],axis = 1)

	desc.columns = ["Description", "Percent_Identity", "Query"]

	assert desc["Query"].is_unique, "Description dataframe contains duplicate queries. Please remove these before continuing."
	
	return desc

####################
# HEADER FUNCTIONS #
####################

def concatenate_reference_strings(header):
	reference = ""
	for row in range(len(header)):
		reference += header.iloc[row,2]
	return reference.upper()

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

#######################
# ALIGNMENT FUNCTIONS #
#######################

def split_query_column_from_alignment(alignment):
	return alignment[0].str.split(n = 1, pat = " +", expand = True)

def create_list_of_groups(query_alignment, desc):
	query_alignment['next_query'] = query_alignment[0].shift(-1)
	query_alignment['last_alignment_row'] = np.where( 
	( ((query_alignment[0] == desc['Query'].iloc[-1]) & (query_alignment['next_query'] == desc['Query'].iloc[0])) | (query_alignment['next_query'].isna()) ),
	True, False)

	df_split_indexes = query_alignment.index[query_alignment['last_alignment_row']]
	list_n_rows = np.diff(df_split_indexes)
	list_n_rows = np.insert(list_n_rows, 0, df_split_indexes[0]+1)

	list_of_groups = np.arange(len(df_split_indexes))
	groups = np.repeat(list_of_groups, list_n_rows)
	return groups

def split_alignment_dataframes_by_group(query_alignment):
	return [d.iloc[:, 0:2].reset_index(drop = True) for _, d in query_alignment.groupby(['groups'])]

def split_columns(df_list, groups, header):
	new_df_list = list()

	reference_length = header.iloc[-1]
	n_nucleotides = 150 #change to snakemake params!!

	for i in range(groups.max() + 1):

		# extract start position column
		df1 = df_list[i][1].str.split(n = 1, pat = "  ", expand = True)

		if i == 0: # if first dataframe, append query id and start position
			new_df_list.append(df_list[i][0])
			new_df_list.append(df1[0])

		# remove extra whitespaces due to differences in start position str length
		col = df1[0].str.len()
		df1['new_col'] = col.max() - col
		x = df1.apply(lambda x: x[1][x["new_col"]:], 1)

		if i < groups.max()	or reference_length % n_nucleotides == 0: # if not last column or last column contains max number of nucleotides
			n = n_nucleotides + 1
		else: # number of nucleotides in last column varies
			n = reference_length % n_nucleotides + 1

		# split remaining columns
		df2 = x.str.split(pat = "", n = n, expand = True)

		if i < groups.max(): # if not last column: don't include end position
			new_df_list.append(df2.iloc[:, 1:-1])
		else: # include end position in final dataframe
			new_df_list.append(df2.iloc[:, 1:])

	return new_df_list


def concatenate_alignment_columns(df_list):
	return pd.concat(df_list, axis=1, join = 'inner')


####################
# CREATE DATAFRAME #
####################

def create_final_dataframe(desc_file, alignment_file, header_file):
	alignment = pd.read_csv(alignment_file, header = None)
	header = create_header(header_file)
	desc = create_description(desc_file)

	query_alignment = split_query_column_from_alignment(alignment)
	groups = create_list_of_groups(query_alignment, desc)
	query_alignment['groups'] = groups
	df_list = split_alignment_dataframes_by_group(query_alignment)
	df_list = split_columns(df_list, groups, header)
	alignment = concatenate_alignment_columns(df_list)
	alignment.columns = header
	df = alignment.merge(desc, how = 'left')

	return df

df = create_final_dataframe(snakemake.input[0], snakemake.input[1], snakemake.input[2])

#######################
# CREATE OUTPUT FILES #
#######################

df.to_csv(snakemake.output[0])

df_nucleotide_changes = df[df["Percent_Identity"] != 100]
df_nucleotide_changes.to_csv(snakemake.output[1])

df_count_nucleotide_changes = df.iloc[:,3:-3].apply(pd.value_counts).fillna(0).filter(regex = '[a-zA-Z]', axis=0).convert_dtypes()
df_count_nucleotide_changes.to_csv(snakemake.output[2])