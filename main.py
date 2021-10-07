from Bio import SeqIO
from itertools import product

START_CODON = 'ATG'
STOP_CODONS = ['TAA', 'TAG', 'TGA']
FILE_FORMAT = 'fasta'
MIN_BP_SIZE = 100
ALL_CODONS = ["TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG", "TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA",
              "TGG", "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", "CAT", "CAC", "CAA", "CAG", "CGT", "CGC",
              "CGA", "CGG", "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG", "AAT", "AAC", "AAA", "AAG", "AGT",
              "AGC", "AGA", "AGG", "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG", "GAT", "GAC", "GAA", "GAG",
              "GGT", "GGC", "GGA", "GGG"]


def split_by_triplets(string_seq):
    triplets = [string_seq[i:i + 3] for i in range(0, len(string_seq), 3)]
    if len(triplets[-1]) < 3:
        triplets.pop()
    return triplets


def split_by_six(string_seq):
    six = [string_seq[i:i+6] for i in range(0, len(string_seq), 6)]
    if len(six) > 0 and len(six[-1]) < 6:
        six.pop()
    return six


def filter_protein_coding_frames_larger_than(frames):
    f = []
    for frame in frames:
        if len(frame) >= MIN_BP_SIZE:
            f.append(frame)
    return f


def get_protein_coding_frames(seq_rec):
    triplets = split_by_triplets(str(seq_rec))
    protein_coding_frames = []
    start_found = 0
    i = 0
    j = 0
    while i < len(triplets):
        if triplets[i] == START_CODON:
            j = i + 1
            while j < len(triplets):
                if triplets[j] in STOP_CODONS:
                    break
                j = j + 1
            protein_coding_frames.append(triplets[i:j + 1])
            i = j + 1
        else:
            i = i + 1
    return protein_coding_frames


def find_stop_start_sequences(str_seq):
    codon_list = split_by_triplets(str_seq)
    i = 0
    longest_frames_between_stop_start = []
    while i < len(codon_list):
        if codon_list[i] in STOP_CODONS:
            start_codons_count = 0
            j = i + 1
            while j < len(codon_list):
                if codon_list[j] == START_CODON:
                    start_codons_count = start_codons_count + 1
                    start_poz = j
                if codon_list[j] in STOP_CODONS:
                    if start_codons_count > 0:
                        longest_frames_between_stop_start.append(codon_list[i:start_poz + 1])
                        i = j - 1
                    break
                j = j + 1
        i = i + 1
    return longest_frames_between_stop_start


def find_codon_freq(str_seq_all_protein_coding_frames_concat):
    codon_list = split_by_triplets(str_seq_all_protein_coding_frames_concat)
    usage_list = []
    for codon in ALL_CODONS:
        usage_list.append((codon, codon_list.count(codon)/len(codon_list)))
    return usage_list


def find_dicodon_freq(all_protein_coding_frames):
    dicodon_list = []
    usage_list = []
    for frame in all_protein_coding_frames:
        dicodons = split_by_six(''.join(frame))
        dicodons.extend(split_by_six(''.join(frame[1:])))
        if len(dicodons) > 0:
            dicodon_list.extend(dicodons)
    all_dicodons = [(i+j) for i in ALL_CODONS for j in ALL_CODONS]
    for dicodon in all_dicodons:
        usage_list.append((dicodon, dicodon_list.count(dicodon)/len(dicodon_list)))
    return usage_list


def get_frames(seq):
    f1 = seq[0:]
    f2 = seq[1:]
    f3 = seq[2:]
    rev_f1 = (seq[::-1])[0:]
    rev_f2 = (seq[::-1])[1:]
    rev_f3 = (seq[::-1])[2:]
    all_frames = [f1, f2, f3, rev_f1, rev_f2, rev_f3]
    return all_frames


def calculate_distances(list_of_tuples_of_file_codon_freq_dicodon_freq):
    phylip_number = len(list_of_tuples_of_file_codon_freq_dicodon_freq)
    result_string = str(phylip_number) + '\n'
    for tuple in list_of_tuples_of_file_codon_freq_dicodon_freq:
        result_string = result_string + tuple[0] + ' '
        for other_tuple in list_of_tuples_of_file_codon_freq_dicodon_freq:
            sum = 0
            for tuple_codon in tuple[1]:
                for other_tuple_codon in other_tuple[1]:
                    if tuple_codon[0] == other_tuple_codon[0]:
                        if tuple_codon[1]/other_tuple_codon[1] >= 5:
                            print(tuple[0] + '   ' + other_tuple[0] + '   ' + str(tuple_codon[0]))
                        sum = sum + abs(tuple_codon[1] - other_tuple_codon[1])
            result_string = result_string + str(sum) + ' '
        result_string = result_string + '\n'
        with open('results_with_codon_usage.txt', 'w') as f:
            f.write(result_string)
            f.write('\n')
    result_string = str(phylip_number) + '\n'
    for tuple in list_of_tuples_of_file_codon_freq_dicodon_freq:
        result_string = result_string + tuple[0] + ' '
        for other_tuple in list_of_tuples_of_file_codon_freq_dicodon_freq:
            sum = 0
            for tuple_dicodon in tuple[2]:
                for other_tuple_dicodon in other_tuple[2]:
                    if tuple_dicodon[0] == other_tuple_dicodon[0]:
                        sum = sum + abs(tuple_dicodon[1] - other_tuple_dicodon[1])
            result_string = result_string + str(sum) + ' '
        result_string = result_string + '\n'
        with open('results_with_dicodon_usage.txt', 'w') as f:
            f.write(result_string)
            f.write('\n')


def util_flatten_list_of_lists(list_of_lists):
    return sum(list_of_lists, [])


def execute(file):
    seq_record = SeqIO.read(file + ".fasta", FILE_FORMAT)
    str_seq = str(seq_record.seq)
    seq_frames = get_frames(str_seq)
    protein_coding_frames = []
    longest_stop_start_frames = []
    filtered_protein_coding_frames = []
    for frame in seq_frames:
        protein_coding_frames.extend(get_protein_coding_frames(frame))
        longest_stop_start_frames.extend(find_stop_start_sequences(frame))
    filtered_protein_coding_frames = filter_protein_coding_frames_larger_than(protein_coding_frames)
    codon_freq = find_codon_freq(''.join(util_flatten_list_of_lists(protein_coding_frames)))
    with open(file + '_res.txt', 'w') as f:
        f.write('PROTEIN CODING FRAMES:\n')
        for pcf in protein_coding_frames:
            f.write(str(pcf) + '\n')
        f.write('\n\n---------------------------------------------------\n\n')
        f.write('LONGEST_STOP_START_FRAMES:\n')
        for lssf in longest_stop_start_frames:
            f.write(str(lssf) + '\n')
        f.write('\n\n---------------------------------------------------\n\n')
        f.write('FILTERED PROTEIN CODING FRAMES LONGER THAN ' + str(MIN_BP_SIZE) + '\n')
        for fpcf in filtered_protein_coding_frames:
            f.write(str(fpcf) + '\n')
    dicodon_freq = find_dicodon_freq(protein_coding_frames)
    p = (file, codon_freq, dicodon_freq)
    return p


if __name__ == '__main__':
    b1 = execute("bacterial1")
    b2 = execute("bacterial2")
    b3 = execute("bacterial3")
    b4 = execute("bacterial4")
    m1 = execute("mamalian1")
    m2 = execute("mamalian2")
    m3 = execute("mamalian3")
    m4 = execute("mamalian4")
    calculate_distances([b1, b2, b3, b4, m1, m2, m3, m4])
    print('COMPLETE')