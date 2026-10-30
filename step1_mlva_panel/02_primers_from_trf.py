#!/usr/bin/env python
import sys, os, getopt 

def usage():
	print("Usage: vntr_from_trf.py -i <input_directory> [-d <delimiter_type>-<value>]")
	print("Example: vntr_from_trf.py -i ../GCF_000262545_DSM_20514 -d inf-6")
	print("  -i, --input: Path to directory containing TRF output files (required)")
	print("  -d, --delimiter: Filter VNTRs by motif size (default: inf-6)")
	print("    inf-X: Keep VNTRs with motif size < X bp")
	print("    sup-X: Keep VNTRs with motif size >= X bp")

try:
	opts, args = getopt.getopt(sys.argv[1:], "i:d:o:", ["input=","delimiter="])
except getopt.GetoptError as err:
	print(f"Error: {err}")
	usage()
	sys.exit(2)

# Initialize variables with defaults
path = None
delimiter = "inf"
delimiter_value = 6

for opt, arg in opts:
	if opt in ("-i", "--input"):
		path = arg
	elif opt in ("-d", "--delimiter"):
		try:
			delimiter = arg.split("-")[0]
			delimiter_value = int(arg.split("-")[1])
		except (IndexError, ValueError):
			print("Error: Delimiter format should be 'inf-VALUE' or 'sup-VALUE'")
			usage()
			sys.exit(2)

# Check if input path was provided
if path is None:
	print("Error: Input path is required!")
	usage()
	sys.exit(2)

# Check if path exists
if not os.path.isdir(path):
	print(f"Error: Directory '{path}' does not exist!")
	sys.exit(2)

files = os.listdir(path)
files = [file for file in files if 'txt' in file]

if not files:
	print(f"Warning: No .txt files found in '{path}'")

dico_comp = {'A':'T','C':'G',"G":"C","T":"A","M":"K","R":"Y","W":"W","S":"S","Y":"R","K":"M","V":"B","H":"D","D":"H","B":"V","X":"X","N":"X",".":".","|":"|"}

def inverComp(seq):
	"""Return the inverse complementary sequence"""
	seq = seq.upper()
	seq_comp = ""
	for nuc in seq:
		seq_comp += dico_comp[nuc]
	seq_comp = "".join(reversed(seq_comp))
	return seq_comp

vntr = []
tmp_left = []
tmp_right = []
left = False
right = False

if path[-1] != '/':
	path = path + '/'

for file in files:
	try:
		with open(path + file) as vntr_file:
			for line in vntr_file.readlines():
				if '<A NAME=' in line:
					loci = line.split('"')[1].replace(",", "_")
					tmp = loci.split("_")
					size_insert = int(tmp[0].split("-")[2]) - int(tmp[0].split("-")[0]) + 1 + 60
					loci = "_".join([tmp[0], str(tmp[3])+"bp", str(size_insert)+"bp", str(int(round(float(tmp[2]))))+""+"U"])

				if left == True:
					tmp_left.append(line)
					if len(line) < 65:
						left = False
						tmp_left = "".join(tmp_left).replace("\r", "").replace("\n", "")[470:495]

				if right == True:
					tmp_right.append(line)
					if len(line) < 65:
						right = False
						tmp_right = inverComp("".join(tmp_right).replace("\r", "").replace("\n", "")[5:30])

				if 'Left flanking sequence' in line:
					left = True
				if 'Right flanking sequence' in line:
					right = True
				if 'Found at' in line and 'loci' in locals():
					vntr.append([loci, tmp_left, tmp_right])
					tmp_left = []
					tmp_right = []
	except IOError as e:
		print(f"Warning: Could not read file {file}: {e}")
		continue

# Apply delimiter filter
if delimiter == "inf":
	vntr = [line for line in vntr if int(line[0].split("_")[1].replace("bp", "")) < delimiter_value]
elif delimiter == "sup":
	vntr = [line for line in vntr if int(line[0].split("_")[1].replace("bp", "")) >= delimiter_value]

# Sort VNTR list by locus start
vntr = [[int(e[0].split("-")[0])] + e for e in vntr]
vntr = [e[1:] for e in sorted(vntr)]

# If there are 2 VNTRs at the same locus, keep the one with the smaller motif
previous = '0--0_0bp_0bp_0U'
vntr2 = []
for v in vntr:
	if v[0].split('-')[0] == previous.split("-")[0] or v[0].split('-')[2].split("_")[0] == previous.split("-")[2].split("_")[0]:
		motif_size_v = int(v[0].split("_")[1].replace('bp', ''))
		motif_size_previous = int(previous.split("_")[1].replace('bp', ''))
		if motif_size_previous > motif_size_v:
			vntr2[-1] = v
	else:
		vntr2.append(v)
	previous = v[0]

# Generate output filename
output_path = f"primers_{os.path.basename(os.path.dirname(path))}_{delimiter}{delimiter_value}bp.txt"

# Write output
with open(output_path, "w") as output:
	for line in vntr2:
		output.write("\t".join(line) + "\n")

print(f"Processed {len(vntr2)} VNTRs (filter: {delimiter}{delimiter_value}bp)")
print(f"Output written to: {output_path}")
