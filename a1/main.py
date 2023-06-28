# Andrew Shen
# BIMM 182 Assignment #1
# 4/5/23 - 4/11/23

# Problem 1
def simple_helloworld():
    """
    Print the phrase "Hello Bioinformatics"
    :return: none
    """
    print("Hello Bioinformatics")


# Problem 3
def cat_func(fasta_in):
    """
    Return string containing header of each sequence and length of sequence
    :param fasta_in: string in fasta format
    :return: string of sequence headers and sequence lengths
    """
    seq_str = ""

    # Iterate through all sequences
    fasta_split = fasta_in.split(">")[1:]
    for seq in fasta_split:
        curr_seq = seq.split("\n")
        seq_str += curr_seq[0] + "\n"

        # Iterate through each line of current sequence
        curr_length = 0
        for i in range(1, len(curr_seq)):
            if len(curr_seq[i]) > 0:
                curr_length += len(curr_seq[i])
        seq_str += str(curr_length) + "\n" + "\n"

    return seq_str


# Problem 4
def filter_func(fasta_in):
    """
    Return string containing only mouse and rat sequences
    :param fasta_in: string in fasta format
    :return: string in fasta formats
    """
    seq_dict = {}

    # Populate dictionary with key: header and value: sequence
    fasta_split = fasta_in.split(">")[1:]
    for seq in fasta_split:
        curr_seq = seq.split("\n")
        seq_dict[curr_seq[0]] = "".join(curr_seq[1:])

    # Filter for only mouse and rat sequences
    str_out = ""
    for k, v in seq_dict.items():
        if "Rattus norvegicus" in k or "Mus musculus" in k or "ratAurA" in k or "PDPK_RAT" in k or "M3KC_RAT" in k:
            # Format back into fasta format, with only 60 characters per line
            str_out += ">" + k + "\n\n"
            curr_index = 0
            while len(v[curr_index:]) > 60:
                str_out += v[curr_index:curr_index+60] + "\n\n"
                curr_index += 60
            str_out += v[curr_index:] + "\n\n"

    return str_out.strip()


# Problem 5
def create_indices(fasta_in):
    """
    Create two index files calculated from the input fasta file
    :param fasta_in: string in fasta format
    :return: None
    """

    # Create file data.seq
    seq_str = ""
    fasta_split = fasta_in.split(">")[1:]
    for seq in fasta_split:
        curr_seq = seq.split("\n")
        seq_str += "".join(curr_seq[1:])
        seq_str += "@"
    seq_str = seq_str[:-1]
    # Write to output file
    seq_out = open("data.seq", "w")
    seq_out.write(seq_str)
    seq_out.close()

    # Create file data.in
    in_str = ""
    curr_index = 0
    fasta_split = fasta_in.split(">")[1:]
    for seq in fasta_split:
        curr_seq = seq.split("\n")
        header_split = curr_seq[0].split("|")
        in_str += header_split[1] + " "
        in_str += str(curr_index) + "\n"
        curr_index += len("".join(curr_seq[1:])) + 1

    in_out = open("data.in", "w")
    in_out.write(in_str)
    in_out.close()


# Problem 6
def get_seq(query, seq_file, in_file):
    """
    Return gi number of database sequence containing query string
    :param query: string input query to be found
    :param seq_file: string format of data.seq index file
    :param in_file: string format of data.in index file
    :return: string gi number
    """
    # Find starting index of query sequence
    start = seq_file.index(query)

    # Format "in_file" into dictionary
    in_list = [val.split() for val in in_file.splitlines()]

    # Use binary search to find what range "start" falls into, inf loop
    low_index = 0
    high_index = len(in_list)
    curr_index = None
    curr_gi = None
    while True:
        # If reached the final value
        if curr_index == (high_index + low_index) // 2:
            break
        else:
            curr_index = (high_index + low_index) // 2
        curr_gi = in_list[curr_index][0]
        curr_loc = in_list[curr_index][1]
        # If value is in lower half
        if start < int(curr_loc):
            high_index = curr_index
        # If value is in higher half
        elif start > int(curr_loc):
            low_index = curr_index
    return curr_gi


# Main function
if __name__ == '__main__':

    # Read in data files
    with open("datafile.txt") as file:
        datafile_in = file.read()
    with open("data.seq") as file:
        dataseq_in = file.read()
    with open("data.in") as file:
        datain_in = file.read()

    # Problem 1
    # simple_helloworld()

    # Problem 3
    cat_out = cat_func(datafile_in)
    #  print(cat_out.strip())

    # Problem 4
    filter_out = filter_func(datafile_in)
    # print(filter_out)

    # Problem 5
    create_indices(datafile_in)

    # Problem 6
    query_test = "MHIQITDFGTAKVLSPDS"
    query_out = get_seq(query_test, dataseq_in, datain_in)
    # print(query_out)
