import subprocess
import tempfile
import os
from matplotlib import pyplot as plt
from Bio import pairwise2
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import motifs
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Align.AlignInfo import SummaryInfo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Align.Applications import ClustalOmegaCommandline
from sklearn.cluster import DBSCAN
from sklearn.manifold import MDS


def two_seq_alignment(sequences, output_file):
    seq1 = Seq("ATCGCTAGCTAGCTGCTAGCTAGCTGCTAGCTAGCTAGCTAGCTAGCTA")
    seq2 = Seq("GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")
    with open(output_file, "w") as file:
        alignments = pairwise2.align.globalms(Seq(sequences[0]), Seq(sequences[1]), 1, -1, -2, -2)
        file.write(pairwise2.format_alignment(*alignments[0]))
        # for alignment in alignments:
        #     print("aaa")
        #     file.write(pairwise2.format_alignment(*alignment))


def alignment_clustalo_muscle(in_file, out_file, algorithm_is_muscle=True):
    if algorithm_is_muscle:
        alignment_with_muscle(in_file, out_file)
    else:
        alignment_with_clustalo(in_file, out_file)


def alignment_with_clustalo(in_file, out_file):  # write in file
    if len(get_sequences_from_fasta_file(in_file)) == 1:
        alignment_with_muscle(in_file, out_file)
    else:
        clustalo_exe = r"clustal-omega-1.2.2-win64\clustalo.exe"
        clustalo_cline = ClustalOmegaCommandline(clustalo_exe, infile=in_file, outfile=out_file, force=True)
        stdout, stderr = clustalo_cline()


def alignment_with_muscle(in_file, out_file):  # write in file
    muscle_exe = r"muscle.exe"
    # muscle commands:  muscle -align test_seq.txt -output output.muscle
    #                   muscle -align test_seq.txt -stratified -output output.muscle
    #                   muscle -disperse output.muscle
    command = [muscle_exe, "-align", in_file, "-output", out_file]  # command constructing
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    # muscle_cline = MuscleCommandline(muscle_exe, input=in_file, out=out_file)
    # print(muscle_cline)
    # stdout, stderr = muscle_cline()


def mds_format(table):
    matrix = [x for x in table]
    mds = MDS(dissimilarity='precomputed')
    transformed_matrix = mds.fit_transform(matrix)
    return transformed_matrix


def clustering(table, e):
    clusters = DBSCAN(eps=e, min_samples=2, metric='precomputed').fit(table)
    groups = dict()
    arr = dict()
    for seq_id, cluster_id in zip(table.names, clusters.labels_):
        arr[seq_id] = cluster_id
        if groups.get(cluster_id) is not None:
            groups[cluster_id].append(seq_id)
        else:
            groups[cluster_id] = [seq_id]
    print("Number of clusters: ", len(groups.keys()))
    print("Clusters: ", end="")
    for key in groups.keys():
        print(*groups[key], sep=", ", end=" || ")
    print()
    # for seq_id, cluster_id in zip(table.names, clusters.labels_):
    #     print(f"{seq_id} belongs to cluster {cluster_id}")
    print("Core cluster is: ", end="")
    for i in clusters.core_sample_indices_:
        print(table.names[i], end=" ")
    print()
    return arr, groups


def get_sequences_from_fasta_file(file):
    sequences = []
    for seq in SeqIO.parse(file, "fasta"):
        sequences.append(seq.seq)
    return sequences


IUPAC_codes = {'A': ['A'],
               'C': ['C'],
               'G': ['G'],
               'T': ['T'],
               'R': ['A', 'G'],
               'Y': ['C', 'T'],
               'W': ['A', 'T'],
               'S': ['C', 'G'],
               'K': ['G', 'T'],
               'M': ['A', 'C'],
               'B': ['C', 'G', 'T'],
               'D': ['A', 'G', 'T'],
               'H': ['A', 'C', 'T'],
               'V': ['A', 'C', 'G'],
               'N': ['A', 'C', 'G', 'T']}


def get_consensus_letter(sequences, index):
    alphabet = dict([('A', 0), ('C', 0), ('G', 0), ('T', 0)])
    for seq_index in range(len(sequences)):
        letter = sequences[seq_index][index]
        if letter == '-':
            continue
        elif IUPAC_codes.get(letter) is not None:
            for let in IUPAC_codes.get(letter):
                alphabet[let] += 1
    sorted_alphabet = sorted(alphabet.items(), key=lambda x: -x[1])
    letters = [sorted_alphabet[0][0]]
    for i in range(1, len(sorted_alphabet)):
        if sorted_alphabet[i][1] == sorted_alphabet[i - 1][1]:
            letters.append(sorted_alphabet[i][0])
        else:
            break
    return list(IUPAC_codes.keys())[list(IUPAC_codes.values()).index(letters)]


def get_consensus_sequence(input_file, output_file, built_in=False):
    if built_in:
        output_file.write(">\n" + str(
            SummaryInfo(AlignIO.read(input_file.name, "fasta")).dumb_consensus(threshold=0.5, ambiguous='N')) + "\n")
    else:
        sequences = get_sequences_from_fasta_file(input_file.name)
        print(sequences)
        consensus_str = ""
        for i in range(len(sequences[0])):
            consensus_str += get_consensus_letter(sequences, i)
        output_file.write(">\n" + consensus_str + "\n")


def get_consensus_from_clusters(clusters, seq_file_name, res_file, algorithm_is_muscle):
    intermediate_res = open("clusters_consensus.txt", "w")  # file with consensus sequences for each cluster
    sequence_dict = SeqIO.to_dict(SeqIO.parse(seq_file_name, "fasta"))
    tmp = tempfile.NamedTemporaryFile(mode="w", delete=False)
    print(sequence_dict)
    for cluster_id in clusters.keys():
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as cluster:
            for sequence_id in clusters.get(cluster_id):
                print(sequence_id)
                cluster.write(">\n" + str(sequence_dict[sequence_id].seq) + "\n")
        print(len(list(SeqIO.parse(cluster.name, "fasta"))))
        alignment_clustalo_muscle(cluster.name, tmp.name, algorithm_is_muscle)
        print(len(list(SeqIO.parse(tmp.name, "fasta"))))
        get_consensus_sequence(tmp, intermediate_res, False)

    # motif = motifs.create(list(SeqIO.parse(tmp.name, "fasta")))
    # intermediate_res.write(">\n" + str(motif.consensus + "\n"))

    intermediate_res.close()
    alignment_clustalo_muscle("clusters_consensus.txt", tmp.name, algorithm_is_muscle)
    get_consensus_sequence(tmp, res_file, False)

    # motif = motifs.create(list(SeqIO.parse(tmp.name, "fasta")))
    # consensus = str(motif.consensus)

    tmp.close()
    os.remove(tmp.name)
    # return consensus


colors = ['r', 'g', 'b', 'c', 'w']


def print_2d_distance_matrix(table, clusters, name):
    matrix = [x for x in table]
    mds = MDS(dissimilarity='precomputed')
    transformed_matrix = mds.fit_transform(matrix)
    # plt.scatter(transformed_matrix[:, 0], transformed_matrix[:, 1])
    plt.title(name + " distance table")
    for i, (x, y) in enumerate(zip(transformed_matrix[:, 0], transformed_matrix[:, 1])):
        c = clusters.get(table.names[i])
        plt.scatter(x, y, color='black' if c == -1 else colors[c])
        # plt.text(x, y, table.names[i], va='bottom', ha='center', fontsize='large')


def main():
    # print(consensus_seq("output.muscle"))
    # alignment_with_clustalo("seq.txt", "output.clustalo")
    # alignment_with_muscle("seq.txt", "output.muscle")
    print("alignment done")
    # muscle_alignment = AlignIO.read("output8001200.muscle", "fasta")
    clustalo_alignment = AlignIO.read("all_sequences800-1200.muscle", "fasta")

    calc = DistanceCalculator('identity')  # distance matrices

    clustalo_distance_table = calc.get_distance(clustalo_alignment)
    # muscle_distance_table = calc.get_distance(muscle_alignment)

    clustalo_clusters = clustering(clustalo_distance_table, 0.2)
    # muscle_clusters = clustering(muscle_distance_table, 0.2)
    # for i in range(1, 10):
    #     print("epsilon: ", i / 10)
    #     clustering(muscle_distance_table, i / 10)
    #     clustering(clustalo_distance_table, i / 10)

    # fig = plt.figure(2, (8, 4))
    # fig.add_subplot(121)  # making images
    print_2d_distance_matrix(clustalo_distance_table, clustalo_clusters[0], "clustalo")
    # fig.add_subplot(122)
    # print_2d_distance_matrix(muscle_distance_table, muscle_clusters[0], "muscle")
    # fig.subplots_adjust(wspace=0.2, hspace=0.5)

    plt.show()

    res = open("consensus8001200.txt", "w")  # file with final consensus sequences
    get_consensus_from_clusters(clustalo_clusters[1], "all_sequences800-1200.txt", res, False)
    # get_consensus_from_clusters(muscle_clusters[1], "seq.txt", res, False)
    res.close()


main()



# two_seq_alignment(get_sequences_from_fasta_file("consensus.txt"), "consensus_pairwise.txt")
# alignment_with_clustalo("consensus.txt", "consensus.clustalo")
# alignment_with_muscle("consensus.txt", "consensus.muscle")


