from Bio import Entrez
import time
startTime = time.time()

mail = input("Your mail address (you will be notified by NCBI in case of request-issues): ")
term = input("Search term (in case of multiple search queries, separate terms by 'or'): ")
search_period = input("Search period in number of days from now (1 month = 30, 1 year = 365): ")

Entrez.email = mail  # Always tell NCBI who you are.

search_handle = Entrez.esearch(db="nucleotide", term=term, usehistory="y", idtype="acc", datetype="pdat", reldate=search_period)  

search_results = Entrez.read(search_handle)
search_handle.close()

count = int(search_results["Count"])
webenv = search_results["WebEnv"]
query_key = search_results["QueryKey"]

print("Total number of records:", count)

answer = input(f"Do you want to continue with the download for search term '{term}' with {count} hits (Y/n)? ").lower()

if answer == "y" or answer == "yes":
    print("Continuing with download... \n")

    batch_size = 500
    out_handle = open(snakemake.output[0], "w")

    for start in range(0, count, batch_size):
        end = min(count, start + batch_size)
        print("Downloading records %i to %i" % (start + 1, end))
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

else:
    print("Terminating...")