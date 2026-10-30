#!/usr/bin/env python
import sys, os, getopt 

#example : vntr_from_trf.py [output directory TRF]

try:		#check if correct args
	opts, args = getopt.getopt(sys.argv[1:], "i:d:o:", ["input=","delimiter="])
except getopt.GetoptError as err:
	usage()
	sys.exit(2)

for opt, arg in opts: #get args given by user
		if opt in ("-i", "--input"):
			path=arg
		elif opt in ("-d","--delimiter") :
			delimiter=arg.split("-")[0]
			delimiter_value=int(arg.split("-")[1])

files = os.listdir(path)
files=[file for file in files if 'txt' in file]

#print files

dico_comp = {'A':'T','C':'G',"G":"C","T":"A","M":"K","R":"Y","W":"W","S":"S","Y":"R","K":"M","V":"B","H":"D","D":"H","B":"V","X":"X","N":"X",".":".","|":"|"}

def inverComp (seq) :						#return the inversed complementary sequence : TTCGA -> TCGAA
	seq = seq.upper()                  		#acgt -> ACGT
	seq_comp = ""
	for nuc in seq :                   		#nuc stand for nucleotide
		seq_comp += dico_comp[nuc]  		#use of the dictionnary to have the complementary nucleotide
	seq_comp = "".join(reversed(seq_comp))	#we reverse the sequence	
	return (seq_comp)

vntr = []
tmp_left=[]
tmp_right=[]
left=False
right=False

if path[-1]!= '/' : path=path+'/'
for file in files : 
	with open(path+file) as vntr_file :
		for line in vntr_file.readlines() :
			if '<A NAME=' in line : 
				loci = line.split('"')[1].replace(",","_")
				tmp=loci.split("_")
				size_insert=int(tmp[0].split("-")[2]) - int(tmp[0].split("-")[0]) + 1 + 60 #on ajoute 1 car la derniere base est inclue, 60 car on prend en compte la taille des primers et le decalage de 5 soit : 2*(25+5)
				loci = "_".join([tmp[0],str(tmp[3])+"bp",str(size_insert)+"bp",str(int(round(float(tmp[2]))))+"U"])
				#print loci

			if left==True :
				tmp_left.append(line)
				if len(line)<65 : 
					left=False
					tmp_left="".join(tmp_left).replace("\r","").replace("\n","")[470:495]
					#print tmp_left

			if right==True :
				tmp_right.append(line)
				if len(line)<65 : 
					right=False
					tmp_right=inverComp("".join(tmp_right).replace("\r","").replace("\n","")[5:30])
					#print tmp_right

			if 'Left flanking sequence' in line : 
				left = True
			if 'Right flanking sequence' in line : 
				right = True
			if 'Found at' in line and 'loci' in locals() :
				vntr.append([loci,tmp_left,tmp_right])
				tmp_left=[]
				tmp_right=[] 

if delimiter =="inf" :
	vntr = [line for line in vntr if int(line[0].split("_")[1].replace("bp","")) < delimiter_value] 
elif delimiter =="sup" :
	vntr = [line for line in vntr if int(line[0].split("_")[1].replace("bp","")) >= delimiter_value]

#sort of vntr list by locus start
vntr=[[int(e[0].split("--")[0])]+e for e in vntr]
vntr=[e[1:] for e in sorted(vntr)]


#if there are 2 vntr at the same locus, we keep the vntr with the smaller motif
previous='0--0_0bp_0bp_0U'
vntr2=[]
for v in vntr : 
	if v[0].split('-')[0]==previous.split("-")[0] or v[0].split('-')[2].split("_")[0]==previous.split("-")[2].split("_")[0]: #if same locus start or locus end
		motif_size_v=int(v[0].split("_")[1].replace('bp',''))
		motif_size_previous=int(previous.split("_")[1].replace('bp',''))
		if motif_size_previous > motif_size_v :
			vntr2[-1]=v
	
	else : vntr2.append(v)
	previous=v[0]

output_path="primers_"+os.path.dirname(path).replace("TRF_output_","")+"_"+delimiter+str(delimiter_value)+"bp.txt"

with open(output_path,"w") as output :
	for line in vntr2 :
		output.write("\t".join(line)+"\n")
"""
for line in vntr2 :
	print "\t".join(line)
"""

