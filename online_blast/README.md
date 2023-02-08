# Online Blast
This workflow uses txt files from online nucleotide BLAST searches to create tables of aligned sequences and their description.

To use this tool, you will need to download a specific format from NCBI after performing the nucleotide BLAST search:  

![How to download a correctly formatted txt file from an NCBI BLAST search](imgs/blastformat.png?raw=true)

1. Go to the Alignments tab  
2. Change alignment view to 'Flat query-anchored with dots for identities'  
3. Change line length to 150  
4. Press the Download button in the right most part of the settings bar  
5. Press the 'Text' option to download a txt file containing descriptions and alignments  
6. Move the downloaded txt file to a folder named 'input' within this 'local_blast' folder  
7. Repeat this process until the input folder contains all the txt files you wish to investigate
8. Run the workflow entering the following command in bash: ```snakemake -j1```

![Required file format](imgs/fileformat.png?raw=true)

In order for this tool to work, the downloaded txt file format must look exactly like the one shown in the figure above.
In case NCBI changes the format of the txt file, this tool will not work.
Specifically this tool requires:  
1. There must be a '%' sign somewhere in the description rows, and nowhere else  
2. Accession no. must be in the last column of the descriptions
3. Per. Ident. must be in the third to last column
4. Any descriptive text must be in first in the description, followed by 10 additional columns
5. The alignment column must be the 3rd column, and must contain at least 4 dots, i.e. 4 nucleotides must be identical to the reference
6. The first column of the alignment header must be named 'Query'
