from Bio import Entrez, SeqIO

Entrez.email = "slavailyin2004@gmail.com"

# Search query
min_len = 1001
max_len = 1200
# search_query = "Klebsiella pneumoniae blaKPC AND (" + str(min_len) + ":" + str(max_len) + "[SLEN])"
search_query = "Homo sapiens KRT19"


start = 0
batch_size = 20

seq_number = 0

with open("Homo sapiens KRT19.txt", "w") as output_handle:
    while True:
        # Use the Entrez.esearch function to search for the gene
        # search_handle = Entrez.esearch(db="nucleotide", term=search_query, retmax=batch_size, retstart=start)
        search_handle = Entrez.esearch(db="gene", term=search_query, retmax=1, retstart=start)
        search_results = Entrez.read(search_handle)
        search_handle.close()

        id_list = search_results["IdList"]

        if not id_list:
            break

        for record_id in id_list:
            fetch_handle = Entrez.efetch(db="nucleotide", id=record_id, rettype="fasta", retmode="text")
            print(record_id)
            record = SeqIO.read(fetch_handle, "fasta")

            try:
                output_handle.write(">" + record.id + " " + record.description + "\n" + str(record.seq) + "\n")
            except Exception as e:
                print(e)
            fetch_handle.close()

        start += batch_size

print("Downloaded " + str(start) + "sequences")





