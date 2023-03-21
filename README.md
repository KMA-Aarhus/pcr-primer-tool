# PCR primer tool
This tool is created for the Department of Clinical Microbiology at Aarhus University Hospital to aid in monitoring of mutations in PCR primer sites. 

There are two different snakemake workflows: One using an online blast on NCBI and another which performs a local blast search using a downloaded database.

The local search option allows detection of more than 5000 hits and is more user friendly. The upper limit on number of hits can easily be updated.
However, if you require a very large database, use the online version. 

### How to run
Clone the repository.

Activate snakemake, decide on whether you are using the local blast or downloading a blast result. Then navigate to the appropriate subdir.

#### online_blast
To run using online_blast, follow the guide here:
  
[online_blast](https://github.com/KMA-Aarhus/pcr-primer-tool/tree/main/online_blast) for at generere påkrævede filer. Naviger dernæst til online_blast og kør:
```
cd online_blast
snakemake -j1 --use-conda --conda-frontend mamba
```
Note that the first time the tool is run, a number of packages will be installed. For this, we recommend using mamba. 

#### local_blast
To run local_blast, insert your query sequence in the "input" directory (We're working on enabling multiple input sequences)

```
cd local_blast
snakemake -j1 --use-conda --conda-frontend mamba
```
Note that the first time the tool is run, a number of packages will be installed. For this, we recommend using mamba. 

You will be asked to enter the following: 
- "Email" (in case of NCBI support).
- "Search term": organism(s) of interest. Here, it is highly recommended to be as precise as possible. For example, instead of "influenza A", write "H1N1 influenza A AND (viruses[filter] AND host='homo sapiens'". This minimises the risk of incorrect hits. Note you can always check your search result on the NCBI website.
- "Search period in number of days from now": The neccesary timespan. 
The tool will then tell you how many records are returned. For a quick sanity check, ask youself the following:
- Is this a reasonable number for my search?
- How many total bases am I expecting? (Expected genome size)x(total number of records). Fx 100.000 genomes of size 10kbp will take up 1gb of space while  5.000 genomes of size 5Mbp will take up 25gb of space. 

If you are happy with the above, type "Yes" at the "Do you want to continue with the download" prompt. The database will be downloaded and the analysis run automatically.

The database is downloaded to database/db.fasta. To BLAST another sequence from the same organism(s), simply change the query.txt.
If you are BLASTing against a different organism, rename or delete db.fasta. We recommend keeping the database for 1 month in case of additional queries.

#### Output
When the tool is done, you can find 3 files in the "output" directory:
- nucleotide_changes_table.csv: Summarises up all positions in the query sequence including which and how many changes have been found in the database.
- nucleotide_changes_table.csv: Summarises all positions in the query sequence as well as all BLAST hits containing changes.
- alignments_with_nucleotide_changes.csv: Summarises all positions in the query sequence plus all BLAST hits.
