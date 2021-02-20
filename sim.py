from Bio import SeqIO
import glob
import pandas as pd


#Get a read count list.
Samples_readcount=[]

#Get file names of sample files
filenames= glob.glob("/home/hema/Desktop/simulated_reads/*.fasta")
#print filenames

# Iterate over sample files
for fn in filenames:
	#Read sample file
	records = list(SeqIO.parse(fn, "fasta"))
	#Store read count of a sample file into a dict
	occurrences_dict = {}
	for seqrecord in records:
		idkeep, rest = seqrecord.id.split('/',1)
		occurrences_dict[rest] = occurrences_dict.get(rest, 0) + 1
		#print occurrences_dict
	#Append dict of read count of sample file into a list
	Samples_readcount.append(occurrences_dict)
#print (Samples_readcount[0])
Finalmatrix = pd.DataFrame(Samples_readcount)
#print Finalmatrix
k = Finalmatrix.T
#print k
k.to_csv('test1.csv')

