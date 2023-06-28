# Andrew Shen
# BIMM 182 Assignment #2
# 4/13/23 - 4/26/23

# Import packages
import random
from matplotlib import pyplot as plt
import numpy as np
import time


# Problem 1
def loc_al(seq_files, match, mismatch, indel, align):
    """
    Perform local alignment and output score and alignment
    :param seq_files: string in fasta format
    :param match: integer match score
    :param mismatch: integer mismatch score
    :param indel: integer indel score
    :param align: boolean to return alignment or not
    :return: integer local alignment score and string alignment
    """

    # Format input files
    seq_split = seq_files.splitlines()
    seq1 = seq_split[0]
    seq2 = seq_split[1]

    # Set scoring values
    match_score = match
    mismatch_score = mismatch
    indel_score = indel

    # Implement local alignment
    max_index = [0, 0]
    dp = [[0] * (len(seq2) + 1) for i in range(len(seq2) + 1)]
    backtrack = [[0] * (len(seq2) + 1) for i in range(len(seq2) + 1)]

    # Calculate DP array
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):

            deletion = dp[i - 1][j] + indel_score
            insertion = dp[i][j - 1] + indel_score
            # curr_diagonal is either match or mismatch
            if seq1[i-1] == seq2[j-1]:
                curr_diagonal = dp[i - 1][j - 1] + match_score
            else:
                curr_diagonal = dp[i - 1][j - 1] + mismatch_score
            zero = 0

            # Populate DP and backtrack arrays
            dp[i][j] = max(deletion, insertion, curr_diagonal, zero)
            # 1=deletion, 2=insertion, 3=match
            if dp[i][j] == deletion:
                backtrack[i][j] = 1
            elif dp[i][j] == insertion:
                backtrack[i][j] = 2
            elif dp[i][j] == curr_diagonal:
                backtrack[i][j] = 3
            # Update max_index
            if dp[i][j] >= dp[max_index[0]][max_index[1]]:
                max_index = [i, j]

    # Perform backtracking
    str1 = ""
    str2 = ""
    curr = max_index
    while True:
        i = curr[0]
        j = curr[1]
        val = backtrack[i][j]
        if val == 1:
            str1 = seq1[i - 1] + str1
            str2 = "-" + str2
            curr = [i - 1, j]
        elif val == 2:
            str1 = "-" + str1
            str2 = seq2[j - 1] + str2
            curr = [i, j - 1]
        elif val == 3 or val == 4:
            str1 = seq1[i - 1] + str1
            str2 = seq2[j - 1] + str2
            curr = [i - 1, j - 1]
        # Ready to return
        if backtrack[curr[0]][curr[1]] == 0:
            # print("Starting index: " + str([curr[0], curr[1]]))
            break

    # Return score and alignment
    str_out = "Score: " + str(dp[max_index[0]][max_index[1]]) + "\n" + "Length: " + str(len(str1))
    if align:
        # Format for blast alignment
        str_out += "\n\n"
        blast_format = ""
        for i in range(len(str1)):
            if str1[i] == str2[i]:
                blast_format += "|"
            else:
                blast_format += " "
        # Format by  60-char chunks
        index = 0
        while len(str1[index:]) > 60:
            str_out += str1[index:index+60] + "\n" + blast_format[index:index+60] + "\n" + str2[index:index+60] + "\n\n"
            index += 60
        str_out += str1[index:] + "\n" + blast_format[index:] + "\n" + str2[index:]
    # Temp comment to do Problem #2
    # return str_out
    return len(str1)


# Problem 2
def random_dna(seq_num, seq_length):
    """
    Generate random DNA sequences of a certain length
    :param seq_num: integer number  of sequences
    :param seq_length: integer length of sequences
    :return: string of random DNA sequences and summary
    """
    # Generate random sequences
    nucleotides = ["A", "C", "G", "T"]
    sequences = []
    nucleotide_freq = {"A": 0, "C": 0, "G": 0, "T": 0}
    for i in range(seq_num):
        curr_seq = ""
        for j in range(seq_length):
            curr_nuc = random.choice(nucleotides)
            nucleotide_freq[curr_nuc] += 1
            curr_seq += curr_nuc
        sequences.append(curr_seq)

    # Format random sequences and summary for output
    str_out = "\n".join(sequences)
    str_out += "\n\nSummary of Observed Nucleotide Frequencies\n"
    for k, v in nucleotide_freq.items():
        str_out += k + ": " + str(v) + "\n"

    return str_out.strip()


def plot_histogram(data):
    """
    Plot a histogram of one set of data
    :param data: integer list of data
    :return: None, printed histogram
    """
    # Convert input data into a numpy array
    data.sort()
    data_np = np.array(data)

    # Plot histogram
    min_val = data_np[0]
    max_val = data_np[-1]
    fig, ax = plt.subplots(figsize=(8, 5))
    bins = [min_val + i*(max_val-min_val)/10 for i in range(10)]
    ax.hist(data_np, bins)
    plt.title("Histogram of Lengths of Local Alignment")
    plt.xlabel("Lengths in bases")
    plt.ylabel("Number of occurrences")
    plt.show()


# Problem 4
def loc_al_linear(seq_files, match, mismatch, indel, align):
    """
    Perform local alignment and output score and alignment
    :param seq_files: string in fasta format
    :param match: integer match score
    :param mismatch: integer mismatch score
    :param indel: integer indel score
    :param align: boolean to return alignment or not
    :return: integer local alignment score and string alignment
    """

    # Format input files
    seq_split = seq_files.splitlines()
    seq1 = seq_split[0]
    seq2 = seq_split[1]

    # Set scoring values
    match_score = match
    mismatch_score = mismatch
    indel_score = indel

    # CALCULATE ending coordinate
    print("Calculating ending coordinate...")
    max_index = [0, 0]
    max_score = 0
    row2 = [0] * (len(seq2) + 1)
    for i in range(1, len(seq1) + 1):
        row1 = row2
        row2 = [0]
        if i%1000 == 0:
            print(i)
        for j in range(1, len(seq2) + 1):
            deletion = row1[j] + indel_score
            insertion = row2[j - 1] + indel_score
            if seq1[i - 1] == seq2[j - 1]:
                curr_diagonal = row1[j - 1] + match_score
            else:
                curr_diagonal = row1[j - 1] + mismatch_score
            zero = 0
            row2.append(max(deletion, insertion, curr_diagonal, zero))
            if row2[-1] >= max_score:
                max_index = [i, j]
                max_score = row2[-1]

    # CALCULATE starting coordinate
    print("Calculating starting coordinate...")
    max_index_reversed = [0, 0]
    max_score_reversed = 0
    seq1_reversed = seq1[::-1]
    seq2_reversed = seq2[::-1]
    row2_reversed = [0] * (len(seq2_reversed) + 1)
    for i in range(1, len(seq1_reversed) + 1):
        row1_reversed = row2_reversed
        row2_reversed = [0]
        if i%1000 == 0:
            print(i)
        for j in range(1, len(seq2_reversed) + 1):
            deletion = row1_reversed[j] + indel_score
            insertion = row2_reversed[j - 1] + indel_score
            if seq1_reversed[i - 1] == seq2_reversed[j - 1]:
                curr_diagonal = row1_reversed[j - 1] + match_score
            else:
                curr_diagonal = row1_reversed[j - 1] + mismatch_score
            zero = 0
            row2_reversed.append(max(deletion, insertion, curr_diagonal, zero))
            if row2_reversed[-1] >= max_score:
                max_index_reversed = [i, j]
                max_score_reversed = row2_reversed[-1]

    # CALCULATE backtracking matrix from start to ending coordinates
    print("Calculating backtracking matrix...")
    start = [len(seq1)-max_index_reversed[0], len(seq1)-max_index_reversed[1]]
    end = max_index

    dp = [[0] * (len(seq2) + 1) for i in range(len(seq2) + 1)]
    backtrack = [[0] * (len(seq2) + 1) for i in range(len(seq2) + 1)]

    # Calculate DP array
    for i in range(start[0], end[0]+1):
        for j in range(start[1], end[1]+1):

            deletion = dp[i - 1][j] + indel_score
            insertion = dp[i][j - 1] + indel_score
            # curr_diagonal is either match or mismatch
            if seq1[i - 1] == seq2[j - 1]:
                curr_diagonal = dp[i - 1][j - 1] + match_score
            else:
                curr_diagonal = dp[i - 1][j - 1] + mismatch_score
            zero = 0

            # Populate DP and backtrack arrays
            dp[i][j] = max(deletion, insertion, curr_diagonal, zero)
            # 1=deletion, 2=insertion, 3=match
            if dp[i][j] == deletion:
                backtrack[i][j] = 1
            elif dp[i][j] == insertion:
                backtrack[i][j] = 2
            elif dp[i][j] == curr_diagonal:
                backtrack[i][j] = 3

    # Perform backtracking
    str1 = ""
    str2 = ""
    curr = end
    while True:
        i = curr[0]
        j = curr[1]
        val = backtrack[i][j]
        if val == 1:
            str1 = seq1[i - 1] + str1
            str2 = "-" + str2
            curr = [i - 1, j]
        elif val == 2:
            str1 = "-" + str1
            str2 = seq2[j - 1] + str2
            curr = [i, j - 1]
        elif val == 3 or val == 4:
            str1 = seq1[i - 1] + str1
            str2 = seq2[j - 1] + str2
            curr = [i - 1, j - 1]
        # Ready to return
        if backtrack[curr[0]][curr[1]] == 0:
            # print("Starting index: " + str([curr[0], curr[1]]))
            break

    # Return score and alignment
    str_out = "Score: " + str(dp[end[0]][end[1]]) + "\n" + "Length: " + str(len(str1))
    if align:
        # Format for blast alignment
        str_out += "\n\n"
        blast_format = ""
        for i in range(len(str1)):
            if str1[i] == str2[i]:
                blast_format += "|"
            else:
                blast_format += " "
        # Format by  60-char chunks
        index = 0
        while len(str1[index:]) > 60:
            str_out += str1[index:index + 60] + "\n" + blast_format[index:index + 60] + "\n" + str2[
                                                                                               index:index + 60] + "\n\n"
            index += 60
        str_out += str1[index:] + "\n" + blast_format[index:] + "\n" + str2[index:]

    # print("Max Score: " + str(max_score))
    # print("Max Score Reversed: " + str(max_score_reversed))
    # print("Start: " + str(start))
    # print("End: " + str(end))

    return str_out


# Main function
if __name__ == '__main__':

    # Read in files
    with open("p1seqs.txt") as file:
        p1_read = file.readlines()
        p1_files = p1_read[1] + p1_read[4]
    with open("p4pairs.txt") as file:
        p4_read = file.readlines()
        p4_files = p4_read[1] + p4_read[4]

    # -----PROBLEM 1-----
    p1_reg = loc_al(p1_files, 1, -10, -1, True)

    # -----PROBLEM 2-----

    # Generate 500 pairs of 1000bp sequences
    pair1 = random_dna(500, 1000)
    pair2 = random_dna(500, 1000)
    pair1_list = pair1.splitlines()[0:500]
    pair2_list = pair2.splitlines()[0:500]

    # Align 500 sequence pairs of length 1000bp
    # Alignment P1
    length_list_p1 = []
    for i in range(len(pair1_list)):
        if i%50 == 0:
            print(i)
        curr_pair = pair1_list[i] + "\n" + pair2_list[i]
        length_list_p1.append(loc_al(curr_pair, 1, -30, 0, False))
    print(length_list_p1)
    print()
    alignment_p1_lengths = ['1354', '1353', '1361', '1354', '1354', '1357', '1347', '1353', '1355', '1356', '1354', '1354', '1347', '1347', '1354', '1348', '1351', '1354', '1356', '1350', '1346', '1351', '1356', '1354', '1345', '1359', '1359', '1357', '1343', '1355', '1355', '1339', '1351', '1359', '1352', '1356', '1365', '1347', '1351', '1349', '1348', '1356', '1361', '1352', '1354', '1352', '1344', '1338', '1361', '1347', '1358', '1357', '1349', '1343', '1343', '1358', '1345', '1358', '1360', '1352', '1355', '1351', '1348', '1353', '1352', '1346', '1350', '1347', '1362', '1351', '1343', '1350', '1349', '1357', '1340', '1349', '1342', '1351', '1354', '1356', '1354', '1358', '1352', '1356', '1350', '1344', '1351', '1354', '1352', '1352', '1350', '1362', '1350', '1345', '1351', '1345', '1353', '1354', '1354', '1356', '1360', '1354', '1344', '1351', '1349', '1354', '1349', '1355', '1349', '1351', '1364', '1355', '1346', '1358', '1357', '1348', '1343', '1348', '1350', '1350', '1357', '1348', '1351', '1351', '1362', '1342', '1347', '1351', '1354', '1352', '1349', '1353', '1348', '1347', '1360', '1345', '1363', '1345', '1351', '1340', '1349', '1351', '1357', '1350', '1352', '1350', '1352', '1355', '1354', '1352', '1351', '1351', '1344', '1364', '1352', '1355', '1352', '1358', '1350', '1354', '1353', '1363', '1355', '1346', '1349', '1355', '1352', '1349', '1356', '1349', '1360', '1344', '1348', '1349', '1337', '1355', '1354', '1356', '1350', '1348', '1356', '1344', '1346', '1355', '1349', '1359', '1350', '1357', '1353', '1350', '1345', '1355', '1358', '1356', '1356', '1355', '1353', '1347', '1356', '1348', '1356', '1339', '1359', '1350', '1347', '1357', '1355', '1355', '1346', '1349', '1350', '1354', '1357', '1351', '1351', '1346', '1352', '1362', '1353', '1352', '1359', '1354', '1357', '1351', '1351', '1353', '1352', '1347', '1352', '1364', '1354', '1354', '1350', '1350', '1346', '1355', '1354', '1348', '1351', '1353', '1349', '1347', '1357', '1348', '1353', '1354', '1350', '1353', '1356', '1343', '1347', '1342', '1344', '1349', '1347', '1357', '1351', '1349', '1343', '1350', '1349', '1354', '1347', '1348', '1356', '1362', '1348', '1351', '1361', '1352', '1353', '1348', '1350', '1347', '1357', '1349', '1357', '1345', '1345', '1359', '1357', '1359', '1355', '1355', '1348', '1353', '1347', '1352', '1346', '1358', '1350', '1351', '1349', '1358', '1359', '1363', '1353', '1348', '1353', '1360', '1350', '1345', '1356', '1353', '1352', '1359', '1349', '1349', '1352', '1347', '1347', '1352', '1352', '1351', '1358', '1348', '1344', '1355', '1353', '1349', '1350', '1355', '1349', '1350', '1342', '1353', '1359', '1357', '1349', '1357', '1348', '1343', '1357', '1349', '1355', '1355', '1346', '1362', '1345', '1349', '1357', '1348', '1345', '1357', '1345', '1358', '1360', '1356', '1352', '1351', '1346', '1350', '1355', '1352', '1355', '1355', '1358', '1356', '1353', '1356', '1349', '1357', '1355', '1362', '1344', '1355', '1362', '1356', '1350', '1356', '1348', '1342', '1360', '1360', '1347', '1353', '1345', '1364', '1351', '1354', '1348', '1351', '1355', '1349', '1346', '1353', '1356', '1350', '1358', '1347', '1345', '1361', '1349', '1353', '1359', '1354', '1358', '1350', '1349', '1355', '1347', '1350', '1352', '1359', '1349', '1352', '1354', '1356', '1355', '1343', '1353', '1358', '1350', '1353', '1357', '1355', '1362', '1352', '1356', '1347', '1349', '1351', '1356', '1349', '1359', '1353', '1354', '1354', '1359', '1350', '1348', '1346', '1346', '1350', '1352', '1350', '1347', '1357', '1359', '1341', '1349', '1352', '1353', '1349', '1352', '1340', '1353', '1351', '1349', '1366', '1350', '1356', '1346', '1352', '1353', '1355', '1352', '1350', '1359', '1358', '1353', '1348', '1347', '1356', '1358', '1357', '1339', '1356', '1342', '1352', '1351', '1361', '1362', '1346', '1348', '1353', '1347', '1342', '1349', '1348', '1350', '1354', '1347', '1360', '1360', '1346', '1357', '1355', '1347', '1352', '1362', '1352', '1346', '1361', '1353', '1353', '1353', '1350', '1351', '1358']
    alignment_p1_lengths_int = [int(v) for v in alignment_p1_lengths]
    print(sum(alignment_p1_lengths_int) / len(alignment_p1_lengths_int))

    # Alignment P2
    length_list_p2 = []
    for i in range(len(pair1_list)):
        if i%50 == 0:
            print(i)
        curr_pair = pair1_list[i] + "\n" + pair2_list[i]
        length_list_p2.append(loc_al(curr_pair, 1, -30, -20, False))
    print(length_list_p2)
    alignment_p2_lengths = ['9', '9', '10', '10', '10', '10', '11', '9', '9', '9', '9', '9', '10', '10', '10', '9', '11', '11', '10', '11', '9', '11', '9', '10', '10', '9', '11', '10', '9', '9', '8', '12', '9', '9', '8', '9', '8', '11', '10', '11', '9', '10', '9', '9', '10', '10', '10', '12', '10', '8', '9', '10', '10', '9', '9', '9', '10', '10', '13', '10', '8', '10', '10', '9', '10', '10', '9', '9', '10', '10', '10', '11', '9', '9', '10', '9', '9', '10', '10', '9', '11', '10', '10', '9', '10', '12', '10', '11', '9', '8', '9', '10', '9', '9', '10', '8', '10', '9', '9', '9', '10', '9', '9', '8', '9', '10', '9', '10', '8', '9', '11', '9', '10', '10', '8', '9', '11', '9', '9', '11', '11', '10', '10', '12', '10', '9', '9', '11', '9', '9', '11', '11', '9', '9', '9', '9', '9', '9', '9', '10', '9', '10', '9', '9', '10', '8', '9', '9', '10', '9', '9', '8', '12', '11', '9', '9', '10', '9', '11', '9', '13', '9', '10', '9', '9', '9', '10', '10', '11', '11', '8', '11', '10', '9', '11', '11', '10', '10', '8', '9', '10', '9', '9', '9', '10', '9', '9', '12', '10', '12', '10', '10', '10', '12', '10', '8', '10', '10', '10', '10', '9', '11', '11', '9', '10', '9', '9', '9', '9', '9', '9', '9', '11', '10', '9', '9', '9', '10', '9', '10', '9', '11', '9', '10', '10', '9', '9', '11', '10', '9', '9', '9', '8', '10', '11', '9', '10', '11', '10', '10', '8', '10', '10', '8', '11', '10', '10', '8', '12', '9', '10', '9', '12', '9', '11', '9', '9', '9', '10', '9', '9', '10', '10', '9', '8', '12', '9', '9', '10', '10', '10', '10', '9', '9', '9', '11', '10', '12', '9', '10', '10', '13', '9', '9', '9', '9', '10', '9', '9', '9', '10', '9', '11', '11', '12', '10', '9', '10', '10', '10', '9', '9', '9', '11', '9', '9', '10', '9', '12', '9', '9', '9', '9', '9', '10', '10', '10', '12', '9', '9', '9', '10', '10', '9', '14', '9', '10', '10', '12', '9', '10', '10', '10', '9', '9', '9', '8', '10', '10', '9', '10', '10', '9', '9', '9', '8', '9', '9', '9', '8', '9', '10', '9', '11', '9', '10', '10', '11', '9', '10', '9', '9', '9', '11', '10', '10', '9', '10', '11', '10', '9', '9', '10', '10', '9', '9', '10', '9', '10', '9', '9', '8', '9', '9', '9', '9', '10', '10', '9', '9', '10', '9', '10', '10', '10', '11', '9', '9', '8', '9', '11', '9', '12', '9', '12', '10', '9', '9', '9', '10', '9', '10', '11', '10', '9', '9', '10', '9', '9', '9', '10', '8', '9', '10', '9', '9', '11', '10', '9', '9', '10', '10', '8', '9', '9', '9', '11', '13', '10', '10', '9', '11', '11', '9', '10', '9', '9', '11', '10', '12', '8', '9', '9', '9', '9', '9', '9', '13', '12', '9', '8', '11', '10', '10', '10', '9', '10', '10', '9', '10', '9', '12', '8', '10', '10', '10', '10', '8', '9', '10', '9', '10', '9', '10', '12', '9', '12', '11', '11', '11', '11', '9', '8', '9', '10', '9', '10', '9', '10', '9']
    alignment_p2_lengths_int = [int(v) for v in alignment_p2_lengths]
    print(sum(alignment_p2_lengths_int) / len(alignment_p2_lengths_int))

    # Plot histograms of local alignment lengths with P1 and P2
    plot_histogram([int(v) for v in alignment_p1_lengths])
    plot_histogram([int(v) for v in alignment_p2_lengths])

    # Confirm hypothesis by plotting multiple random pairs
    # For P1
    start_time = time.time()

    print("FOR P1...")
    for val in [50, 100, 250, 500, 1000]:
        print("Running lengths for 100 sequences of length %s..." % val)
        pair1 = random_dna(100, val).splitlines()[0:100]
        pair2 = random_dna(100, val).splitlines()[0:100]
        lengths = []
        for i in range(len(pair1)):
            curr_pair = pair1[i] + "\n" + pair2[i]
            lengths.append(loc_al(curr_pair, 1, -30, 0, False))
        print(lengths)
        print(sum(lengths)/len(lengths))
        print()

    print()
    print("----------------")
    print()

    # For P2
    print("FOR P2...")
    for val in [50, 100, 250, 500, 1000]:
        print("Running lengths for 100 sequences of length %s..." % val)
        pair1 = random_dna(100, val).splitlines()[0:100]
        pair2 = random_dna(100, val).splitlines()[0:100]
        lengths = []
        for i in range(len(pair1)):
            curr_pair = pair1[i] + "\n" + pair2[i]
            lengths.append(loc_al(curr_pair, 1, -30, -20, False))
        print(lengths)
        print(sum(lengths) / len(lengths))
        print()

    end_time = time.time()
    print(end_time - start_time)

    # -----PROBLEM 3-----
    print("FOR Problem 3...")
    avg_lengths = []
    for val in [-30, -20, -10, -7.5, -5, -2.5, -2, -1.5, -1, -0.5, -0.33, -0.25, 0]:
        print("Running lengths for 1000 sequences of length 250 for val of %s..." % val)
        pair1 = random_dna(1000, 100).splitlines()[0:1000]
        pair2 = random_dna(1000, 100).splitlines()[0:1000]
        lengths = []
        for i in range(len(pair1)):
            curr_pair = pair1[i] + "\n" + pair2[i]
            lengths.append(loc_al(curr_pair, 1, val, val, False))
        # print(lengths)
        # print(sum(lengths)/len(lengths))
        avg_lengths.append(sum(lengths)/len(lengths))
        # print()
    print(avg_lengths)

    # Plot the points
    x_vals = np.array([-30, -20, -10, -7.5, -5, -2.5, -2, -1.5, -1, -0.5, -0.33, -0.25, 0])
    y_vals = np.array([6.289, 6.231, 6.256, 6.301, 6.353, 7.623, 10.38, 20.018, 82.129, 109.205, 113.616, 114.746, 136.342])
    plt.plot(x_vals, y_vals, "o")
    plt.title("Mean lengths of optimal alignment for varying indel/mismatch penalties")
    plt.xlabel("Indel/Mismatch Penalty")
    plt.ylabel("Mean Lengths of Optimal Alignment")
    plt.show()

    # -----PROBLEM 4-----
    start_time = time.time()
    p1_out_linear = loc_al_linear(p1_files, 1, -10, -1, True)
    end_time = time.time()
    print(p1_out_linear)
    print(end_time - start_time)
    # got score=56 with max_index=[21228,3648] and 546 seconds, max_index_reversed=[21104,3531]
    # Start: [21104,3531]
    # End: [21228,3648]
    # Size of backtracking matrix: 124 x 117
    start_time = time.time()
    p4_out = loc_al_linear(p4_files, 1, -10, -1, True)
    print("----------------")
    print(p4_out)
    print("----------------")
    end_time = time.time()
    print(end_time - start_time)