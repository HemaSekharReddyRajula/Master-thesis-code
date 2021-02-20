import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
import scipy.stats as stats
import sys
from math import sqrt

# opening a file and read it 
file = open('FILE_SAMPLE1_MAP.txt', 'r')
foo =file.readlines()
#print foo

#spliting values and save it in a list a.
a=[]
for values in foo:
	a.append(values.split('\t'))
#print a
del a[0]

#Dictionary in which key as patient id and value as filename and respective barcode.
b={}
for i in range(0,len(a)):
	s=a[i][1].split('-')[2]
	if s not in b:
		b[s]=[]
		for j in range(0,len(a)):
			if s in a[j][1]:
				filename = 'rsem.genes.normalized_results'
				if filename in a[j][0]:
					#print filename
					b[s].append((a[j][0],a[j][1]))
#print b

#Dictionary with matched samples.
i = 0
c={}
for key, val in b.iteritems():
	
	if len(val) == 2:
		#print ({key: val})
		
		c.update({key: val})
		i = i + 1
c.pop('A9CF')
#print c
#print i
#print len(c)
#print "c['A9CF']: ", c['A9CF']
#print c.keys()
#print c.values()

#Opening a RNAseq file, read and store in a list a1.
file1 = open('Liver1.txt', 'r')
foo1 =file1.readlines()
#print foo1

a1=[]
for values1 in foo1:
	a1.append(values1)
#print a1
#print a1[0]

#Matching the filename to the RNAseq data(which has same filename with path).
e=[]
for i in a1:
	for j in c.values():
		for n in j:
			for p in n:
				#print p
			
				if p in i:
			#print i
					i.strip()
					#print i.strip()
					e.append(i.strip())
#print e
#print len(e)

#Dictionary consists of key as filename, barcode and value as geneid, normalized count. 

g={}
for j in c.values():
	for n in j:
		k = n
		for va in e:
			for w in k:
				if w in va:
			
					file2 = open( va , 'r')
					foo2 = file2.readlines()
					#print foo2
					
					p = []
					for values in foo2:
						values.strip()
						p.append(values.strip().split('\t'))
				
					#print p
					#break
					
			
					#g.update({k: v})
					g[k]=p
					#break
#print g.keys()
#print len(g.keys())
#print g.values()[:10]

#keys are moving inside to the respective list of values, keeping it in single list.
h=[]

for i in g.keys():
	h1=[]
	h1.append(list(i))
	for i in g[i]:
	    h1.append(i)
	h.append(h1)

#print h[0][:10]

#print len(h[i])

#print len(h[0])
#h[0].pop(1)
#print (h[0][:10])


#Header1 has literally gene_id and (patient id-01A/11A) in a sequential order. 
header1=[]
header1.append('gene_id')

h1=[]
h2=[]
h3=[]
for j in h:
	if j[0][1].split('-')[3] == '01A':
		h1.append(j[0][1].split('-')[2] +'-'+ j[0][1].split('-')[3])
	elif j[0][1].split('-')[3] == '11A':
		h2.append(j[0][1].split('-')[2] +'-'+ j[0][1].split('-')[3])
	elif j[0][1].split('-')[3] == '02A':
		h3.append(j[0][1].split('-')[2] +'-'+ j[0][1].split('-')[3])
		
h1.extend(h2)
h1.extend(h3)
header1.extend(h1)

#print header1
#print len(header1)
#header2 has the literally gene_id and patient id(01A/11A) and at the end '\n'.
f = ', '.join(header1)
#print f, 'cdef'
header2 = ''.join((f,'\n'))

#print header2, 'abcd'
#print len(header2)


#Contents has genes.
contents = []
for i in h[0][2:]:
	con=[]
	con.append(i[0])
	contents.append(con)
#print contents[20:35]
#print len(contents)

#Removing the letters after |.
conten=[]
for i in contents:
	for j in i:
		conten.append(j.split('|'))
#print conten[20:35]

content=[]
for j in conten:
	con11=[]
	con11.append(j[0])
	content.append(con11)
#print content[20:35]
#print len(content)
#contents1 has normalized count values(100 files) with respect to header1 patient id.
contents1=[]
for i in header1[1:]:
	for m in h:
		if i in m[0][1]:
			
			con1=[]
			for j in m[2:]:
				con1.append(j[1])
			contents1.append(con1)
#print contents1[99][:10]
#print (header1[100])
#print header1
#print len(header1)
#print len(contents1[0])
#print h[0][:10]
#print len(contents)
#print len(contents1)
#print contents1[0][:10]

#y has matched genes and normalized count values.
y=[]
for i in range(0,len(content)):
	k=content[i]
	#print k
	for j in range(0, len(contents1)):
		l=contents1[j][i]
		k.append(l)
	y.append(k)
#print y[40]
#print len(y)
#del y[0:29]

#Removing the genes which has '?'.
y2=[]
for i in y:
	if i[0] != '?':
		y2.append(i)
#print(y2[0])
#print len(y2)
#print type(y2[0])


#Separation of geneids and values.
y3 = []
y4 = []
for i in y2:
	y3.append(i[0])
	y4.append(i[1:])

#converting string into float.

n = []
for item in y4:
	
	new=[]
	for item1 in item:
		new.append(float(item1))
	n.append(new)
#print n[3]


n4a=np.mean(n[3])
#print n4a

na = len(n[3])
#print na
mean = sum(n[3]) / na
#print mean
sd = sqrt(sum((x-mean)**2 for x in n[3]) / na)
#print sd
#print np.mean(n[3])
#print np.std(n[3])

n1 = []
n2 = []
for item in n:
	n1.append(item[0:50])
	n2.append(item[50:100])
#print n1[20501:]

n3=[]
for i in range(0, len(n1)):
	n3.append((n1[i], n2[i])) 

#print n3[0]

#Genewise testing.
nn = []
for item in n3:
	new=[]
	for item1 in item:
		new1=[]
		for item2 in item1:
			new1.append(float(item2))
		new.append(new1)
	nn.append(new)
#print nn[20501:]
#print nn


#scipy.stats.mstats.normaltest gives k2 and P values.
nn2=[]
for i in nn:
	nn1=[]
	for j in i:
		array = np.asarray(j)
		norm = stats.normaltest(array)
		nn1.append(norm)
	nn2.append(nn1)
#print nn2
#If the p value is greater than 0.05 then it does pass the normality test.

nn4=[]
for i in nn2:
	nn3=[]
	for j in i:
		if j[1] > 0.05:
			nn3.append(1)	
		else:
			nn3.append(0)
	nn4.append(nn3)
#print nn4

count = 0
count1 = 0
count2 = 0
count3 = 0
for i in nn4:
	if i[0]==0 and i[1]==0:
		count = count +1
	elif i[0]==0 and i[1]==1:
		count1 = count1 +1
	elif i[0]==1 and i[1]==0:
		count2 = count2 +1
	elif i[0]==1 and i[1]==1:
		count3 = count3 +1
#print count
#print count1
#print count2
#print count3
'''
#Histograms of (0,0), (0,1) and  (1,1).

ba1 = n2[20501:]
ba2 = n1[20501:]
ba3 = n2[20500:20501]
ba4 = n1[20500:20501]
ba5 = n2[20492:20493]
ba6 = n1[20492:20493]
#bins = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]

plt.hist( ba1, bins=30, rwidth=0.8)
#plt.hist( ba2, bins=30, rwidth=0.8)
#plt.hist( ba3, bins=30, rwidth=0.8)
#plt.hist( ba4, bins=30, rwidth=0.8)
#plt.hist( ba5, bins=30, rwidth=0.8)
#plt.hist( ba6, bins=30, rwidth=0.8)
plt.xlabel('x')
plt.ylabel('y')
plt.show()


#Samplewise testing.
nsample = []
for i in range(len(n[0])):
	nsample1=[]
	for j in range(len(n)):
		nsample1.append(n[j][i])
	
	nsample.append(nsample1)
#print nsample[0][0:100]
#print len(nsample)

list2=[]
for i in nsample:
	list1=[]
	for j in i:
		
		if j < 1000:
			list1.append(j)
			#print type(j)	
			
		
	list2.append(list1)
#print len(list1)
#print list2[23]
#print len(list2[23])


nn2=[]
for i in nsample:
	
	array = np.asarray(i)
	norm = stats.normaltest(array)
	nn2.append(norm)
#print nn2

nn4=[]
for i in nn2:
	if i[1] > 0.05:
		nn4.append(1)	
	else:
		nn4.append(0)
#print nn4


#Histograms

ba7 = list1[0]
#print len(ba7)
#print max(ba7)
#print sum(ba7)/len(ba7)t( ba7, bins=30, rwidth=0.8

plt.hist( ba7, bins=30, rwidth=0.8)
plt.xlabel('x')
plt.ylabel('y')
plt.show()


#Mean of tumor and normal.
n4=[]
for i in n1:
	mean=sum(i)/len(i)
	n4.append(mean)
#print len(n4)
#print max(n4)
#print min(n4)

n5=[]
for i in n2:
	mean=sum(i)/len(i)
	n5.append(mean)
#print n5[0:10]
#print len(n5)


#Subtracting the values of n4 and n5.
n6=[]
for i in range(0, len(n4)):
	k=n4[i]
	k1=n5[i]
	k2=k-k1
	k3=abs(k2)t( ba7, bins=30, rwidth=0.8
	n6.append(k3)
#print n6
#print len(n6)



#n6a has the greater mean values of tumor than normal.
#n6a1 has the greater values of normal than tumor.
n6a=[]
n6a1=[]
for i in range(0, len(n4)):
	ka=n4[i]
	k1a=n5[i]
	if ka>k1a:
		n6a.append(ka)
	if ka<k1a:
		n6a1.append(k1a)
#print len(n6a)
#print len(n6a1)
#print n6a[0]



#Barchart for n6a(greater mean values of tumor than normal) and n6a1(greater values of normal than tumor).
ba =['Tumor', 'Normal']
y_pos = np.arange(len(ba))
ba1=[len(n6a),len(n6a1)]
#print ba1
plt.bar(y_pos, ba1)
plt.show()


#Histograms of mean of tumor, normal and difference of them.
ba = n4
ba1 = n5
ba2 = n6
bins = [0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000]

#plt.hist( ba, bins, rwidth=0.8)
#plt.hist( ba1, bins, rwidth=0.8)
plt.hist( ba2, bins, rwidth=0.8)
plt.xlabel('x')
plt.ylabel('y')
plt.show()


#n7 contains the values, gene id and s1 has sorted values.
n7=[]
for i in range(0,len(n6)):
	l2=[]
	l=n6[i]
	l1=y2[i][0]
	l2.append(l)
	l2.append(l1)
	n7.append(l2)
#print n7[0]
#print len(n7)
s1=sorted(n7)
#print s1[20482:]

n8=[]
for i in s1[20499:]:
	gene=i[1]
	#print gene
	#gene1=[]
	for j in range(0, len(y2)):
		if gene in y2[j]:
			n8.append(y2[j])
			#print gene1
	#n8.append(gene1)
#print n8
#print len(n8)

#n9 has gene ids and n10 has values of it.
n9=[]
n10=[]
for i in n8:
	n9.append(i[0])
	n10.append(i[1:])
#print n9

nn9=[]
for i in n9:
	nn9.append(int(i))
#print type(nn9[0])

n11=[]
for i in n10:
	n12=[]
	for j in i:
		n12.append(float(j))
	n11.append(n12)
#print n11

n13=[]
n14=[]
for i in n11:
	n13.append(i[0:50])
	n14.append(i[50:100])t( ba7, bins=30, rwidth=0.8
#print len(n13)

n15=[]
for i in range(0, len(n13)):
	n15.append(n13[i])
	n15.append(n14[i])
#print type(n15)
#print len(n15)

#Boxplot of 19 genes which are much difference between tumor and normal expression values.
plt.boxplot(n15)
plt.title('Tumor&Normal')
plt.ylabel('normalized_count')
plt.xlabel(n9)
plt.show()

#contentstring2 has '\n' among 100 files.
contentstring=[]
for i in y2:
	f1 = ', '.join(i)
	contentstring.append(f1)
#print len(contentstring)
#print type(f1)

contentstring2 = '\n'.join(contentstring) + '\n'
#print contentstring2, 'abcd'

#headcontent has header2(gene_id and patient ids) and contentstrings2(gens+100 normalised count values).
headcontent = header2 + contentstring2
#print len(headcontent)
#print (headcontent)


#Writing a file.
fileout = open('out115.csv', 'w')
fileout.write(headcontent)

#Heatmap of 19 genes which are much difference between tumor and normal expression values.
del header1[0]
#print header1

data1 = np.array(n11)
#print len(data1)
ymax= len(data1)

plt.pcolor(data1, cmap=plt.cm.Reds)
plt.ylim(ymax=ymax)
plt.xlabel('Barcode')
plt.ylabel(n9)
plt.show()


import scipy
import pylab
import scipy.cluster.hierarchy as sch
#Heat map with dendrogram of 19 genes which are much difference between tumor and normal expression values.
data2=np.array(n11)
#print type(data2)t( ba7, bins=30, rwidth=0.8
#print len(data2)
#print len(data2[0])

fig = pylab.figure(figsize=(19,100))
ax1 = fig.add_axes([0.09,0.1,0.2,0.6])
fig = pylab.figure(figsize=(19,100))
ax1 = fig.add_axes([0.09,0.1,0.2,0.6])
Y = sch.linkage(data2, method='centroid')
#print type(Y)
#print len(Y)
#print len(Y[0])

Z1 = sch.dendrogram(Y, orientation='right')
#print type(Z1)print len(Z1)
#print len(Z1[0])
#print Z1

axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
#print type(axmatrix)

idx1 = Z1['leaves']
#print type(idx1)
#print len(idx1)

data2 = data2[idx1,:]
#print type(data2)
#print len(data2)
#print len(data2[0])

im = axmatrix.matshow(data2,aspect='auto', cmap=plt.cm.Reds)
axcolor = fig.add_axes([0.91,0.1,0.02,0.6])t( ba7, bins=30, rwidth=0.8
pylab.colorbar(im, cax=axcolor)
fig.show()
fig.savefig('dendrogram1.png', dpi = 300)


#Write a adjusted p values to file and finding a DE and Non DE genes.
pp=[]
for idid in range(3):
	#Runing edgeR method
	f_edge = open('outfile'+str(idid)+'.txt', 'r')
	fo_edge = f_edge.readlines()
	print fo_edge

	a_edge=[]
	for i in fo_edge:
		a_edge.append(i.split())
	#print a_edge

	#Convert list of lists into single list.
	p_edge = reduce(lambda x,y: x+y,a_edge)
	#print p_edge

	#Convert string into float
	mylisttwo = map(float, p_edge)
	#print mylisttwo

	# Limiting float values to two decimals.
	floatlist = []
	for item in mylisttwo:
		floatlist.append("{0:.4f}".format(item))
	#print floatlist

	dictionary = dict(zip(y3, floatlist))
	#print dictionary
	#Converting the values from string to float.
	dictionay1={}
	for key, value in dictionary.iteritems():
		key1 = key
		value1 = float(value)
		dictionay1.update({key1: value1})
	#print dictionay1
t( ba7, bins=30, rwidth=0.8
	#Values less than 0.05
	dict1={}
	count=0
	for ke, valu in dictionay1.iteritems():
		if valu < 0.05:
			dict1.update({ke: valu})
			count = count + 1
	#print dict1
	#print count

	p1=[]
	for key, value in sorted(dict1.iteritems(), key=lambda (k,v): (v,k)):
		 p = "%s: %s" % (key, value)
		 p1.append(p)
	#print p1[0:10]
	pp.append(p1)
#print pp
k=set(pp[0]).intersection(*pp)
#print k
#print len(k)



f=open('bodymap.txt', 'r')
f1=f.readlines()
#print f1

n66=[]
for values in f1:
	values.strip()
	n66.append(values.strip().split('\t'))
#print n66

n65 = []
for item in n66:
	new=[]
	for item1 in item:
		new.append(float(item1))
	n65.append(new)
#print n65

#Histograms of mean of tumor, normal and difference of them.
#ba = n4
#ba1 = n5
ba2 = n6
bins = [0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000]

#plt.hist( ba, bins, rwidth=0.8)
#plt.hist( ba1, bins, rwidth=0.8)
plt.hist( ba2, bins, rwidth=0.8)
plt.xlabel('x')
plt.ylabel('y')
plt.show()


file = open('test_out0.txt', 'r')
file1 = file.readlines()

#print file1

l=[]
for i in file1:
	temp1=[]
	temp=i.strip().split('"')
	#print temp
	temp1.append(temp[1])
	temp1.append(temp[3])
	l.append(temp1)
	#sys.exit()
#print l
l1={}
for i in l:
	l1[i[0]]=i[1]
#print l1.values()


#Genewise testing.
nn = []
for item in l1.values():
	nn.append(float(item))
#print nn
#print len(nn)


#scipy.stats.mstats.normaltest gives k2 and P values.


array = np.asarray(nn)
norm = stats.normaltest(array)
#nn2.append(norm)
#print norm


#If the p value is greater than 0.05 then it does pass the normality test.

nn4=[]
for i in nn2:
	nn3=[]
	for j in i:
		if j[1] > 0.05:
			nn3.append(1)	
		else:
			nn3.append(0)
	nn4.append(nn3)
#print nn4

count = 0
count1 = 0
count2 = 0
count3 = 0
for i in nn4:
	if i[0]==0 and i[1]==0:
		count = count +1
	elif i[0]==0 and i[1]==1:
		count1 = count1 +1
	elif i[0]==1 and i[1]==0:
		count2 = count2 +1
	elif i[0]==1 and i[1]==1:
		count3 = count3 +1
#print count
#print count1
#print count2
#print count3

#Histograms of (0,0), (0,1) and  (1,1).

ba1 = n2[20501:]
ba2 = n1[20501:]
ba3 = n2[20500:20501]
ba4 = n1[20500:20501]
ba5 = n2[20492:20493]
ba6 = n1[20492:20493]
#bins = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]

plt.hist( ba1, bins=30, rwidth=0.8)
#plt.hist( ba2, bins=30, rwidth=0.8)
#plt.hist( ba3, bins=30, rwidth=0.8)
#plt.hist( ba4, bins=30, rwidth=0.8)
#plt.hist( ba5, bins=30, rwidth=0.8)
#plt.hist( ba6, bins=30, rwidth=0.8)
plt.xlabel('x')
plt.ylabel('y')
plt.show()
'''

from math import sqrt
import math

yy3=[]
for i in range(1,1001):
	yy3.append(int(i))
#print yy3

#Write a adjusted p values to file and finding a DE and Non DE genes.
pp=[]
pp1=[]
for idid in range(1,5,1):
	#Runing edgeR method
	f_edge = open('/home/hema/Desktop/matrix_output'+str(idid)+'.txt', 'r')
	fo_edge = f_edge.readlines()
	#print fo_edge

	l=[]
	for i in fo_edge:
		temp1=[]
		temp=i.strip().split('"')
		#print temp
		temp1.append(temp[1])
		temp1.append(temp[3])
		l.append(temp1)
		#sys.exit()
	#print l[0:2000]
	l1={}
	for i in l:
		l1[str(i[0])]=float(i[1])
	#print l1.values()[0:2000]
	#print l1.keys()[0:2000]
	#print l1
	l2=sorted(l1.iteritems(), key=lambda l1:l1[0], reverse=False)
	#print l2
	#print l2[0:2005]
	#print len(l2)
	
	#Convert list of lists into single list.
	#p_edge = reduce(lambda x,y: x+y,l1.values())
	#print p_edge

	#Convert string into float
	mylisttwo = map(float, l1.values())
	#print mylisttwo
	

	
	my=[]
	for i in mylisttwo:
		my1 = "{:f}".format(float(i))
		my.append(my1)
	#print my	
	

	
	# Limiting float values to two decimals.
	floatlist = []
	for item in my:
		floatlist.append("{0:.15}".format(item))
	#print floatlist
	
	dictionary = dict(zip(yy3, my))
	
	l3=[]
	for i in l2:
		if i[1]<0.05:
			l3.append(i[0])
	TP=0
	FP=0	
	for i in l3:
		if i<11:
			TP=float(TP+1)
		else:
			FP=float(FP+1)
	TN=float(900-FP)
	FN=float(100-TP)
	'''
	TP,FP,TN,FN,n=0,0,0,0,0
	for i in l2:
		n=n+1
		if(i[1]<0.05):
			if(n<101):
				TP=float(TP+1)
			else:
				FP=float(FP+1)
		else:
			if(n<101):
				FN=float(FN+1)
			else:
				TN=float(TN+1)
	'''	
	P=float(TP+FN)
	N=float(FP+TN)
	#print TP,FP,TN,FN
	sensitivity=round(TP/P,3)
	#print "Sensitivity:"+str(sensitivity)
	specificity=round(TN/N,3)
	#print "Specificity:"+str(specificity)
	precision=round(TP/((TP+FP)+1),3)
	#print "Precision:"+str(precision)
	recall=round(TP/(TP+FN),3)
	#print "Recall:"+str(recall)
	accuracy=round((TP+TN)/(P+N),3)
	#print "Accuracy:"+str(accuracy)
	f1_score=round((2*TP)/(2*TP+FP+FN),3)
	#print "F1-score:"+str(f1_score)
	mcc=round(((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)+1),3)
	#print "MCC:"+str(mcc)
	#print (str(mcc)+str(recall)+str(f1_score))
	
	negative_predictivevalue=round(TN/((TN+FN)+1),3)
	#print "NPV: "+str(negative_predictivevalue)
	false_positiverate=round(FP/N,3)
	#print "FPR: "+str(false_positiverate)
	false_discoveryrate=round(FP/((FP+TP)+1),3)
	#print "FDR: "+str(false_discoveryrate)
	false_negativerate=round(FN/(FN+TP),3)
	#print "FNR: "+str(false_negativerate)
	#print recall+f1_score+mcc+accuracy+specificity+precision+sensitivity+negative_predictivevalue+false_positiverate+false_discoveryrate+false_negativerate
	#print recall, f1_score, mcc, accuracy, specificity, precision, sensitivity, negative_predictivevalue, false_positiverate, false_discoveryrate, false_negativerate
		
	#print l3
	#print len(l3)
	#print dictionary
	#Converting the values from string to float.
	dictionay1={}
	for key, value in dictionary.iteritems():
		key1 = key
		value1 = float(value)
		dictionay1.update({key1: value1})
	#print dictionay1

	#Values less than 0.05
	dict1={}
	count=0
	for ke, valu in dictionay1.iteritems():
		if valu < 0.05:
			dict1.update({ke: valu})
			count = count + 1
	#print dict1
	#print len(dict1)
	print count
	pp1.append(dict1)

	p1=[]
	for key, value in sorted(dict1.iteritems(), key=lambda (k,v): (v,k)):
		 p = "%s: %s" % (key, value)
		 p1.append(p)
	#print p1[0:10]
	pp.append(p1)
#print len(pp1)
#print pp1[0]
#print len(pp1[1])
#print pp1[0].values()
k=set(pp[0]).intersection(*pp)
#print k
#print len(k)
'''
pp2=[]
count = 0
for i,j in pp1[1].iteritems():
	pp3=[]
	if i < 2000:
		pp3.append((i,j))
	pp2.append(pp3)
	count = count + 1
	#print pp3	
#print count

#for key in sorted(pp1[1].iterkeys()):
#	print "%s: %s" % (key, pp1[1][key])

#Histograms.

ba1 = pp1[0].values()
ba2 = pp1[1].values()
ba3 = pp1[2].values()


#plt.hist( ba1, bins=30, rwidth=0.8)
#plt.hist( ba2, bins=30, rwidth=0.8)
plt.hist( ba3, bins=30, rwidth=0.8)

plt.xlabel('x')
plt.ylabel('y')
plt.show()
'''
