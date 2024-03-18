from Bio import SeqIO
input_file_name = "Homo_sapience_KRT19.txt"
output_file_name = "CDNA_KRT19.txt"

with open(output_file_name, "w") as output_file:
    for record in SeqIO.parse(input_file_name, "fasta"):
        cDNA = record.seq.reverse_complement()
        mRNA = cDNA.reverse_complement().transcribe()
        if "mRNA" in record.description:
            output_file.write(">" + record.description + "\n")
            output_file.write(str(mRNA) + "\n")

#splicing - The process known as splicing removes introns and joins exons.