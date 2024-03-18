from Bio import Entrez, SeqIO

Entrez.email = "slavailyin2004@gmail.com"


class Genome:
    def __init__(self, genome_id, start, stop, minus_strand, exon_index, seq=None, description=None):
        self.genome_id = genome_id
        self.start = start
        self.stop = stop
        self.minus_strand = minus_strand
        self.exon_table_index = exon_index
        self.seq = seq
        self.description = description

    def download_genome(self):
        if self.genome_id is not None:
            fetch_handle = Entrez.efetch(db="nucleotide", id=self.genome_id, rettype="fasta", retmode="text",
                                         seq_start=self.start,
                                         seq_stop=self.stop)
            print(self.genome_id)
            record = SeqIO.read(fetch_handle, "fasta")
            self.description = record.description
            self.seq = record.seq
            fetch_handle.close()
        return self


def get_similar_gene_ids(search_query):
    start = 0
    batch_size = 20
    gene_ids = []
    while True:
        search_handle = Entrez.esearch(db="gene", term=search_query, retmax=batch_size, retstart=start)
        search_results = Entrez.read(search_handle)
        search_handle.close()

        id_list = search_results["IdList"]

        if not id_list:
            break

        gene_ids.extend(id_list)
        start += batch_size
    print("Downloaded " + str(len(gene_ids)) + "gene ids")
    return gene_ids


def get_gene_table(gene_id):
    gene_handle = Entrez.efetch(db="gene", id=gene_id, rettype="gene_table", retmode="text")
    gene_table = gene_handle.read()
    gene_handle.close()
    return gene_table


def parse_genome_table(split_table):
    genome_id = None
    start, stop = int(), int()
    genome_id_not_found = True
    start_exon_table_index = 0
    minus_strand = False
    for line_index, line in enumerate(split_table):
        if genome_id_not_found and line.startswith("Reference"):
            parts = line.split()
            if len(parts) >= 9:
                genome_id = parts[4]
                pos = 0
                for i in range(4, len(parts)):
                    if parts[i] == "from:":
                        pos = i
                        break
                if parts[5] == "(minus" and parts[6] == "strand)":
                    minus_strand = True
                start, stop = int(parts[pos + 1]), int(parts[pos + 3])
                genome_id_not_found = False
        elif line.startswith("Exon table for  mRNA"):
            start_exon_table_index = line_index
            break
    return Genome(genome_id, start, stop, minus_strand, start_exon_table_index)


def download_exons(split_table, genome, exons_file):
    exon_index = genome.exon_table_index
    if genome.genome_id is None and exon_index != 0:
        return
    genome.download_genome()
    print(split_table[exon_index])
    mRNA_name = ""
    try:
        mRNA_name = split_table[exon_index].split()[4]
    except Exception as e:
        print(str(exon_index))
        print(e)
        return
    print(mRNA_name)
    exon_index += 3
    curr_line = split_table[exon_index]
    while curr_line.strip():
        print(exon_index - 3 - genome.exon_table_index)
        exon_boundaries = curr_line.split()[0]
        exon_start, exon_stop = map(int, exon_boundaries.split('-'))
        if genome.minus_strand:
            exon = genome.seq[genome.start - exon_start:genome.start - exon_stop + 1]
        else:
            exon = genome.seq[exon_start - genome.start:exon_stop - genome.start + 1]
        exons_file.write(">" + genome.genome_id + ":" + exon_boundaries + " " +
                         str(exon_index - 3 - genome.exon_table_index) + " " + mRNA_name + "\n")
        exons_file.write(str(exon) + "\n")
        exon_index += 1
        curr_line = split_table[exon_index]


def download_genes(file_name):
    gene_ids = get_similar_gene_ids("Homo sapiens KRT19")
    with open(file_name, "w") as output_handle:
        for gene_id, cnt in zip(gene_ids, range(1)):
            gene_table = get_gene_table(gene_id)
            split_table = gene_table.split('\n')
            print(gene_table)

            genome = parse_genome_table(split_table)
            print(genome.genome_id)
            genome.download_genome()
            print(str(genome.seq))
            # download_exons(split_table, genome, file_name)
            # if genome_id is not None:
            #     fetch_handle = Entrez.efetch(db="nucleotide", id=genome_id, rettype="fasta", retmode="text",
            #                                  seq_start=genome_inf[1],
            #                                  seq_stop=genome_inf[2])
            #     print(genome_id)
            #     record = SeqIO.read(fetch_handle, "fasta")
            #     try:
            #         output_handle.write(">" + record.id + " " + record.description + "\n" + str(record.seq) + "\n")
            #     except Exception as e:
            #         print(e)
            #     fetch_handle.close()


KRT19_gene_id = 3880
with open("Homo_sapiens_KRT19_exons.txt", "w") as output:
    init_table = get_gene_table(KRT19_gene_id)
    print(init_table)
    table = init_table.split('\n')
    curr_genome = parse_genome_table(table)
    print(curr_genome.genome_id)
    download_exons(table, curr_genome, output)
