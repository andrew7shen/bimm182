# Andrew Shen
# BIMM 182 Assignment #3
# 5/14/23 - 5/17/23

# Problem 4
def aho_corasick_construction(patterns):
    """
    Return adjacency list for trie constructed from patterns
    :param patterns: list of pattern strings
    :return: string format for adjacency list
    """
    # Initialize trie as adjacency list with key: start node, value: [[end node, label]]
    trie = {0: []}
    node_index = 1

    # Iterate through all patterns
    for j in range(len(patterns)):
        pattern = patterns[j]
        curr_node = 0
        # Iterate through each char in pattern
        for i in range(len(pattern)):
            curr_symbol = pattern[i]
            # Check for outgoing edges with label curr_symbol
            is_edge = False
            for edge in trie[curr_node]:
                end_node = edge[0]
                label = edge[1]
                if label == curr_symbol:
                    curr_node = end_node
                    is_edge = True
            if not is_edge:
                # Add end label if end of word
                if i == len(pattern)-1:
                    trie[node_index] = [[j, "*"]]
                else:
                    trie[node_index] = []
                trie[curr_node].append([node_index, curr_symbol])
                curr_node = node_index
                node_index += 1
    # Format trie for output
    # trie_str = ""
    # for k, v in trie.items():
    #     if len(v) > 0:
    #         for edge in v:
    #             trie_str += str(k) + " " + str(edge[0]) + " " + edge[1] + "\n"
    # Currently returning the adjacency list for trie, can return trie_str for correct format
    return trie


def aho_corasick_search(text, patterns, trie):
    """
    Return all starting positions where a pattern from "patterns" appears in "text"
    :param text: string text
    :param patterns: list of pattern strings
    :param trie: adjacency list format for trie
    :return: string format for starting positions dictionary
    """
    # Create output dictionary
    count_dict = {}
    for pattern in patterns:
        count_dict[pattern] = 0

    # Run search algorithm
    l = 0
    v = 0
    c = 0
    while c < len(text):
        curr_char = text[c]
        can_transition = False
        next_node = None

        # Check if can transition nodes
        for edge in trie[v]:  # edge example is [1, "P"]
            curr_next = edge[0]
            curr_dest = edge[1]
            if curr_dest == curr_char: # can transition!
                can_transition = True
                next_node = curr_next
                break

        # If success
        if can_transition:
            v = next_node
            c = c + 1
            if trie[v][0][1] == "*":  # reached word end
                pattern_index = trie[v][0][0]
                count_dict[patterns[pattern_index]] += 1
        # If fail
        else:
            c = l + 1
            l = c
            v = 0

    # Convert output dictionary to output string
    str_out = ""
    for k, v in count_dict.items():
        str_out += k + ": " + str(v) + "\n"
    return str_out


if __name__ == '__main__':

    print("Main method running...")

    # -----PROBLEM 2-----
    L = [5, 11, 15, 20, 25, 30, 35, 40]
    m = 10000000
    n = 10000000
    r = 100
    prob2_out = "L\tSpeed-up\tSensitivity\n"
    # Calculate speed-up and sensitivity for all these values
    for val in L:
        speed_up = (m*n)/(m+n+pow(r, 2)*pow(0.25, val)*m*n)
        sensitivity = 1-pow((1-pow(0.85, val)), r-val+1)
        prob2_out += str(val) + "\t" + str(speed_up) + "\t" + str(sensitivity) + "\n"
    print(prob2_out)

    # -----PROBLEM 4-----
    # Dictionary 1
    print("\nReading in files for dictionary1...")
    with open("DNA.txt") as file:  # Length of DNA: 13,652,895
        dna_read = [line.strip() for line in file.readlines()]
        dna = "".join(dna_read[1:])
    with open("queries.txt") as file:
        queries = [val.strip() for val in file.readlines()]
    print("\nStart trie1 construction...")
    trie1 = aho_corasick_construction(queries)
    print("\nStart trie1 search...")
    trie1_output = aho_corasick_search(dna, queries, trie1)

    # Calculate e-values for first dictionary
    keyword_matches = [val.split(": ") for val in trie1_output.splitlines()]
    keyword_out = ""
    n = 13652895
    for val in keyword_matches:
        k = len(val[0])
        e_value = pow(0.25, k)*(n-k+1)
        comp = ""
        if e_value > float(val[1]):
            comp = "less than expectation"
        elif e_value == float(val[1]):
            comp = "identical to expectation"
        elif e_value < float(val[1]):
            comp = "more than expectation"
        keyword_out += val[0] + ": " + val[1] + ", " + str(e_value) + ", " + comp + "\n"
    print(keyword_out)

    # Dictionary 2
    print("\nReading in files for dictionary2...")
    with open("queries2.txt") as file:
        queries2 = [val.strip() for val in file.readlines()]
    print("\nStart trie2 construction...")
    trie2 = aho_corasick_construction(queries2)
    print("\nStart trie2 search...")
    trie2_output = aho_corasick_search(dna, queries2, trie2)

    # Calculate e-values for first dictionary
    keyword_matches = [val.split(": ") for val in trie2_output.splitlines()]
    keyword_out = ""
    n = 13652895
    for val in keyword_matches:
        k = len(val[0])
        e_value = pow(0.25, k) * (n - k + 1)
        comp = ""
        if e_value > float(val[1]):
            comp = "less than expectation"
        elif e_value == float(val[1]):
            comp = "identical to expectation"
        elif e_value < float(val[1]):
            comp = "more than expectation"
        keyword_out += val[0] + ": " + val[1] + ", " + str(e_value) + ", " + comp + "\n"
    print(keyword_out)
