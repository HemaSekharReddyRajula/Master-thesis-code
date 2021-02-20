

import os

for i in range(1,21,1):
	if i<10:
		m = '0'+str(i)
	else:
		m = str(i)
	f_edge = open('/home/hema/Desktop/Methods/simulated_reads_0/sample_'+m+'_1.fasta', 'r')
	f_edge1 = open('/home/hema/Desktop/Methods/simulated_reads_0/sample_'+m+'_2.fasta', 'r')
	os.system('bwa.kit/bwa index /home/hema/Desktop/Methods/simulated_reads_0/chr22.fa')
	os.system('bwa.kit/bwa mem chr22.fa sample_'+m+'_1.fasta sample_'+m+'_2.fasta > aln-pe_sample_'+m+'.sam')
	os.system('bwa.kit/samtools view -bT chr22.fa aln-pe_sample_'+m+'.sam > aln-pe_sample_'+m+'.bam')
	os.system('bwa.kit/samtools sort aln-pe_sample_'+m+'.bam aln-pe_sample_'+m+'_sorted')
	os.system('bwa.kit/samtools index aln-pe_sample_'+m+'_sorted.bam')
	os.system('bwa.kit/samtools idxstats aln-pe_sample_'+m+'_sorted.bam')
	os.system('bwa.kit/samtools idxstats aln-pe_sample_'+m+'_sorted.bam > /home/hema/Desktop/Methods/simulated_reads_0/output/sample_'+m+'-mapped-readcount.txt')



import csv
from collections import OrderedDict

foo=[]

for i in range(1,21,1):
	if i<10:
		m = '0'+str(i)
	else:
		m = str(i)
	file = open('/home/hema/Desktop/Methods/simulated_reads_0/output/sample_'+m+'-mapped-readcount.txt', 'r')
	foo1 = file.readlines()
	foo.append(foo1)
	 
#print foo
#print len(foo)

a=[]
for values in foo:
	a1=[]
	for value in values:
		a1.append(value.split('\t'))
	a.append(a1)
#print a
#print len(a)
#print len(a[0])

b=OrderedDict()

d=[]
for i in range(0,len(a)):
	for j in range(0,20):
		key=a[i][j][0]
		b[key]=a[i][j][2]

	d.append(b.values())
#print b
#print b.keys()
c=[]
for i in b.keys():
	c.append(i)
#print c
f=[]
f.append(c)
e = f + d
#print e


#Writing a file.

w = csv.writer(open("/home/hema/Desktop/Methods/simulated_reads_0/matrix.csv","w"))
w.writerows(zip(*e))
'''
import csv
with open('/home/hema/Desktop/Methods/simulated_reads_0/matrix.csv') as f:
    r = csv.reader(f)
    data = [line for line in r]
with open('/home/hema/Desktop/Methods/simulated_reads_0/matrix0.csv','w') as f:
    w = csv.writer(f)
    w.writerow(['Transcript','S1','S2','S3','S4','S5','S6','S7','S8','S9','S10','S11','S12','S13','S14','S15','S16','S17','S18','S19','S20'])
    w.writerows(data)
'''
