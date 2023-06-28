
# Command line script for local alignment

import sys

# Read in values
seq_files_in = sys.argv[1]
with open(seq_files_in) as file:
    seq_files_read = file.readlines()
    seq_files = seq_files_read[1] + seq_files_read[4]
match = int(sys.argv[3])
mismatch = int(sys.argv[5])
indel = int(sys.argv[7])
if len(sys.argv) == 9 and sys.argv[8] == "-a":
    align = True
else:
    align = False

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
        break

# Return score and alignment
str_out = "Score: " + str(dp[max_index[0]][max_index[1]]) + "\n" + "Length: " + str(len(str1)) + "\n"
if align:
    # Format for blast alignment
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

print(str_out)
