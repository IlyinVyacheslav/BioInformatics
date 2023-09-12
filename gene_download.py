from Bio import Entrez, SeqIO

Entrez.email = "slavailyin2004@gmail.com"

search_query = "Homo sapiens KRT19"

start = 0
batch_size = 20

seq_number = 0

gene_ids = []

while True:
    # Use the Entrez.esearch function to search for the gene
    # search_handle = Entrez.esearch(db="nucleotide", term=search_query, retmax=batch_size, retstart=start)
    search_handle = Entrez.esearch(db="gene", term=search_query, retmax=batch_size, retstart=start)
    search_results = Entrez.read(search_handle)
    search_handle.close()

    id_list = search_results["IdList"]

    if not id_list:
        break

    gene_ids.extend(id_list)
    start += batch_size

print(gene_ids)

with open("Homo sapiens KRT19.txt", "w") as output_handle:
    for gene_id in gene_ids:
        gene_handle = Entrez.efetch(db="gene", id=gene_id, rettype="gene_table", retmode="text")
        gene_table = gene_handle.read()
        gene_handle.close()
        gene_table_lines = gene_table.split('\n')
    
        genome_id = None
        for line in gene_table_lines:
            if line.startswith("Reference"):
                parts = line.split()
                print(gene_id, line)
                if len(parts) >= 9:
                    genome_id = parts[4]
                    pos = 0
                    for i in range(4, len(parts)):
                        if parts[i] == "from:":
                            pos = i
                            break
                    coordinates = int(parts[pos + 1]), int(parts[pos + 3])
                    break
                    
        if genome_id is not None:
            fetch_handle = Entrez.efetch(db="nucleotide", id=genome_id, rettype="fasta", retmode="text",
                                         seq_start=coordinates[0],
                                         seq_stop=coordinates[1])
            print(genome_id)
            record = SeqIO.read(fetch_handle, "fasta")
            try:
                output_handle.write(">" + record.id + " " + record.description + "\n" + str(record.seq) + "\n")
            except Exception as e:
                print(e)
            fetch_handle.close()

print("Downloaded " + str(start) + "sequences")
