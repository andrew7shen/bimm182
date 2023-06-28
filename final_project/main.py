# Andrew Shen
# BIMM182 Final Project
# 5/25/23

def find_slivers(file_path):
    """
    Return list of cycle strings
    :param file_path: path to file with cycle information, ex: "EFM192A_BREAST_amplicon1_cycles.txt"
    :return: list of which paths are cycles
    """

    # Read in file
    with open(file_path) as file_in:
        cycle_info = [line.strip().split("\t") for line in file_in.readlines()]
    file_path_split = file_path.split("/")[1].split("_")
    sample_name = file_path_split[0] + "_" + file_path_split[1]
    amplicon = file_path_split[2]

    # Find all cycles
    segment_rows = []
    cycle_rows = []
    cycle_list = []  # List of all cycles
    for row in cycle_info:
        curr_row = row[0]
        if len(row) == 1 and curr_row != "List of cycle segments":
            curr_row = curr_row.replace('"', "")
            cycle_rows.append(curr_row.split(";"))
        elif len(row) == 5 and row[0] == "Segment":
            segment_rows.append(row)
    for i in range(len(cycle_rows)):
        curr_row = cycle_rows[i]
        curr_cycle_info = curr_row[2][9:]
        if curr_cycle_info[0] != "0":
            cycle_list.append(curr_row)

    # Calculate which segments are slivers
    sliver_ids = []
    sliver_coords = {}
    for segment in segment_rows:
        id = segment[1]
        start = segment[3]
        end = segment[4]
        if int(end)-int(start) < 1000:
            sliver_ids.append(id)
            sliver_coords[id] = [start, end, str(int(end)-int(start))]

    # Check if each cycle contains slivers
    all_slivers = []
    for cycle in cycle_list:
        # Make list of segments in each cycle
        segment_list = [id[:-1] for id in cycle[2][9:].split(",")]
        sliver_segments = []  # Segments that are slivers in each cycle
        for segment in segment_list:
            if segment in sliver_ids and segment not in sliver_segments:
                sliver_segments.append(segment)
        all_slivers.append(sliver_segments)

    # Format all the info in a string
    str_out = ""
    for i in range(len(cycle_list)):
        cycle = cycle_list[i]
        cycle_id = cycle[0]
        cycle_copycount = cycle[1]
        cycle_segments = cycle[2]
        slivers = ",".join(all_slivers[i])
        sliver_locations = ""
        sliver_sizes = ""
        for sliver in all_slivers[i]:
            curr_coords = sliver_coords[sliver]
            sliver_locations += "[" + curr_coords[0] + "," + curr_coords[1] + "]"
            sliver_locations += ","
            sliver_sizes += curr_coords[2] + ","
        str_out += sample_name + "\t" + amplicon + "\t" + "ecDNA" + "\t" + cycle_id + "\t" + cycle_copycount + "\t" + cycle_segments + "\t" + slivers + "\t" + sliver_locations[:-1] + "\t" + sliver_sizes[:-1] + "\n"

    return str_out


if __name__ == '__main__':

    print("Starting main method...")

    # Read in files
    with open("aggregated_results_trimmed.tsv") as file:
        results = [line.strip().split("\t") for line in file.readlines()]
        results[0].insert(0, "Sample index")

    # Find locations of sliver segments in cycles
    final_str = "sample_name\tamplicon\tclassification\tcycle_id\tcopy_count\tsegments\tsliver_segments\tsliver_locations\tsliver_sizes\n"
    all_paths = ["ecDNA_cycle_info/EFM192A_BREAST_amplicon1_cycles.txt", "ecDNA_cycle_info/EFM192A_BREAST_amplicon3_cycles.txt", "ecDNA_cycle_info/ES2_OVARY_amplicon1_cycles.txt", "ecDNA_cycle_info/ES2_OVARY_amplicon3_cycles.txt", "ecDNA_cycle_info/HARA_LUNG_amplicon2_cycles.txt", "ecDNA_cycle_info/HARA_LUNG_amplicon4_cycles.txt", "ecDNA_cycle_info/HARA_LUNG_amplicon10_cycles.txt"]
    for path in all_paths:
        curr_str = find_slivers(path)
        final_str += curr_str

    # Write to output file
    file_out = open("ecDNa_cycles_with_slivers.tsv", "w")
    file_out.write(final_str.strip())
    file_out.close()
    print("Successfully wrote to output file!")
