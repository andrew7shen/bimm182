with open("slivers.csv", "r") as file:
#with open("sample_types.csv", "r") as file:
	file_read = file.read()
	file_in = [val.strip() for val in file_read.split("'")]

count_dict = {}
for sample in file_in:
	if len(sample) > 1:
		sample_type = sample.split("_")[0] + "_" + sample.split("_")[1]
		if sample_type not in count_dict.keys():
			count_dict[sample_type] = 1
		else:
			count_dict[sample_type] += 1

sorted_dict = dict(sorted(count_dict.items(), key=lambda item: item[1]))

for k,v in sorted_dict.items():
	print(k + " " + str(v))