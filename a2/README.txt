Andrew Shen
BIMM 182 Assignment #2
4/13/23 - 4/26/23

List of all files in code submission
1) README.txt: README document that includes details on relevant code/text files and how to run code
2) main.py: primary Python script with functions and code that performs calculations and formatting for Problems 1, 2, 3, and 4 in the homework
3) loc_al.py: Python script that performs local alignment on a pair of sequences through command line
4) p1seqs.txt: input text file of two sequences for Problem 1
5) p4pairs.txt: input text file of two sequences for Problem 4

How to run the code:
For Problem 1:
- Can call loc_al.py with parameters "<seq_files> -m <match> -s <mismatch> -d <indel> -a"
- Include the "-a" tag if full alignment is desired
- loc_al.py command line script or loc_al function in main.py both output "score" and "length" of optimal alignment
- Use "p1seqs.txt" as input for <seq_files>
- Test case use -m of 1, -s of -10, and -d of -1

For Problem 2:
- Can call function random_dna to generate random DNA sequences
- Specify the desired <number of sequences> and <length of sequences> for the generated sequences
- Can call function plot_histogram to plot histogram of an input list of data values
- Can reference main method in main.py script to see usage of these two functions to generate 500 pairs of 1000bp sequences
- Test case use <number of sequences> of 500 and <length of sequences> of 1000
- Can reference main method in main.py script to see computations of l_p1(n) and l_p2(n) with different parameters

For Problem 3:
- Can reference main method in main.py script to see computations for local alignment of varying indel/mismatch values
- Can reference the same main method to see plotting of the average lengths of local alignment

For Problem 4:
- Can call function loc_al_linear to perform local alignment using linear space
- Can reference main method in main.py script to see call of function loc_al_linear to perform local alignment on large paired sequences of length 25,000
- Uses linear space with validation that length of alignment "L" is much smaller than lengths of input sequences
