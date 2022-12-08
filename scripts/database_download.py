from Bio import Entrez
import time
startTime = time.time()

with open('ncbi_mail.txt') as f:
    mail = f.readline()

Entrez.email = mail  # Always tell NCBI who you are.

search_handle = Entrez.esearch(db="nucleotide", term=snakemake.params[0], usehistory="y", idtype="acc", datetype="pdat", reldate=30)  

search_results = Entrez.read(search_handle)
search_handle.close()

count = int(search_results["Count"])
webenv = search_results["WebEnv"]
query_key = search_results["QueryKey"]

print("Total number of records:", count)

batch_size = 100
out_handle = open(snakemake.output[0], "w")

for start in range(0, count, batch_size):
    end = min(count, start + batch_size)
    print("Going to download record %i to %i" % (start + 1, end))
    fetch_handle = Entrez.efetch(
        db="nucleotide",
        rettype="fasta",
        retmode="text",
        retstart=start,
        retmax=batch_size,
        webenv=webenv,
        query_key=query_key,
        idtype="acc",
    )
    data = fetch_handle.read()
    fetch_handle.close()
    out_handle.write(data)
out_handle.close()

executionTime = round(time.time() - startTime)
print("Execution time in seconds:", executionTime)