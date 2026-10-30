#!/usr/bin/env python3
# =============================================================================
# correction_motif_vntr.py  --  PYTHON 3 PORT, NOT VERIFIED AGAINST THE ORIGINAL
#
# The published results were produced by correction_motif_vntr.py in this same
# folder, which is Python 2 and will not run under the Python 3.10 environment
# the rest of this pipeline specifies. No Python 2 interpreter was available to
# diff the two, so this port is mechanical and UNVERIFIED. If you depend on it,
# check its output against the Python 2 original before trusting it.
#
# Exactly these changes were made, and nothing else:
#   1. fractions.gcd -> math.gcd        (removed from fractions in Python 3.9)
#   2. reduce imported from functools   (a builtin in Python 2)
#   3. print statements -> print()
#   4. three zip() results wrapped in list(): they are iterated more than once,
#      and zip returns an iterator in Python 3
#   5. int((insert_size-60)/pattern_found) -> //, so negative numerators floor
#      as they did in Python 2 rather than truncating toward zero
#   6. mlva_output.next() -> next(mlva_output)
#   7. csv.writer opened with newline="" to stop blank rows on some platforms
# =============================================================================
import sys, os, csv
from math import gcd
from functools import reduce

#ce script a pour but de verifier que la taille du motif est correct
#pour cela on compare la difference de taille d'insert entre les souches ainsi que le nombre de motif retrouve et la taille du motif
#on utilise pour cela le pgcd (plus grand commun diviseur) pour verifier

file_handle=open(sys.argv[1],"r")
file = csv.reader(file_handle,delimiter=",")

output_path=sys.argv[1].split('.')[0]+"_corrected.csv"

dico_vntr={}
dico_modified={}

first_columns=[line[:2] for line in file]
file_handle.seek(0)
file = [line[2:] for line in file] #on retire 'key' et 'access'
file=list(zip(*file))
for i,line in enumerate(file) :
	dico_vntr[line[0]]=set([int(e) for e in line[1:] if e!=''])

def distance(x,y) :
	return abs(x-y)

#return a set of distance between the values inside a list (values of the list are compared against each other)
def motif_size_search(values_list) :
	res= set()
	for i,value in enumerate(values_list) :
		res |= set([distance(value,e) for e in [size for size in values_list if size!=value]])
	if len(res)==0 : res=values_list
	return res 

retained=[]
modified=[]
modified2=[]
modified3=[]
discarded=[]
for key in dico_vntr.keys() :
	insert_size=int(key.split('_')[2].replace('bp',''))
	pattern_found = reduce(gcd, motif_size_search(dico_vntr[key])) #gcd= greatest common divisor
	sizeU=int((insert_size-60)//pattern_found) #-60 because we remove the size of the paired end primers
	if int(key.split('_')[1].replace("bp",""))!= pattern_found :

		modified_primer="_".join([key.split('_')[0],str(pattern_found)+"bp",key.split('_')[2],str(sizeU)+"U"])

		if pattern_found != 1 and pattern_found in motif_size_search(dico_vntr[key]) and pattern_found!=int(insert_size):
			#print insert_size, pattern_found
			#print key, sorted(dico_vntr[key])
			#print key, sorted(motif_size_search(dico_vntr[key]))
			#print pattern_found,"\n"
			
			modified.append((key,modified_primer))

		elif pattern_found != 1 : 
			#sorted(dico_vntr[key])
			#print key, sorted(dico_vntr[key])
			#print key, sorted(motif_size_search(dico_vntr[key]))
			#print pattern_found,"\n"
			modified2.append((key,modified_primer))
			

		elif pattern_found==1 :

			tmp=[]
			for size in sorted(dico_vntr[key]) :
				if size-1 not in tmp : tmp.append(size)  
			found_reduced=tmp 

			#print key, sorted(dico_vntr[key])
			#print key, found_reduced
			#print key, sorted(motif_size_search(found_reduced))
			#print reduce(gcd, motif_size_search(found_reduced)),"\n"
			modified3.append((key,"_".join([key.split('_')[0],str(pattern_found)+"bp",key.split('_')[2],str(sizeU)+"U"])))

		else:
			#print key, sorted(dico_vntr[key])
			#print key, sorted(motif_size_search(dico_vntr[key]))
			#print pattern_found,"\n"
			discarded.append(key)

	else :
		#print key, sorted(dico_vntr[key])
		#print key, sorted(motif_size_search(dico_vntr[key]))
		#print pattern_found,"\n"
		retained.append(key)

print('validated primers:')
for e in retained :
	print(e)
print('\n')
print('modified primers:')
for e in modified :
	print(e)
print('\n')
print('primers sans distance egal au pgcd:')
for e in modified2 :
	print(e)
print('\n')
print('primers avec pgcd egal a 1:')
for e in modified3 :
	print(e)
print('\n')
print('discarded primers')
for e in discarded :
	print(e)
file_handle.seek(0)
file=list(zip(*first_columns))+[line for line in file if line[0] in retained] #+[line for line in file if line[0] in modified]
file=list(zip(*file))

out = csv.writer(open(output_path,"w",newline=""), delimiter=',',quoting=csv.QUOTE_ALL)
for line in file :
	out.writerow(line)

new_small=[]

with open('loci_retained.txt',"w") as loci_retained :
	for locus in retained :
		loci_retained.write(locus+"\n")
	for locus in modified :
		if int(locus[1].split("_")[1].replace("bp","")) >= 12 :
			loci_retained.write(locus[1]+"\n")
		else :
			new_small.append(locus)

with open('loci_modified.txt',"w") as loci_modified :
	for locus in modified :
		loci_modified.write("\t".join(locus)+"\n") 

file_handle=open(os.path.basename(sys.argv[1]).replace("predicted_PCR_size_table_","").replace("_filtered","_output"),'r')
mlva_output=csv.reader(file_handle,delimiter=",")

dico_reads={}
next(mlva_output)
for line in mlva_output :
	dico_reads[line[1]]=(line[8],line[10])

with open('primers_corrected.txt',"w") as primers_corrected :
	for locus in retained :
		primers_corrected.write("\t".join([locus]+list(dico_reads[locus]))+"\n")
	for locus in modified :
		if int(locus[1].split("_")[1].replace("bp","")) >= 12 :
			primers_corrected.write("\t".join([locus[1]]+list(dico_reads[locus[0]]))+"\n")

with open('primers_new_small.txt',"w") as primers_small :
	for locus in new_small :
		primers_small.write("\t".join([locus[1]]+list(dico_reads[locus[0]]))+"\n")