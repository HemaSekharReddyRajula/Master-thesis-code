########## Question:4 ############
d={}
for i in range(1, 9):
	d[i] = i*i
#print d

########## Question:5 ############
values=12,34
l=list(values)
t=tuple(values)
#print l
#print t

########## Question:6 ############
'''
class InputOutString(object):
    def __init__(self):
        self.s = ""

    def getString(self):
        self.s = raw_input()

    def printString(self):
        #print self.s.upper()

strObj = InputOutString()
strObj.getString()
strObj.printString()
'''
########### Question:7 ##########
import math
c=50
h=30
items= 100,150,180
value=[]
for d in items:
	value.append(str(int(round(math.sqrt(2*c*float(d)/h)))))
#print ','.join(value)


def putNumbers(n):
    i = 0
    while i<n:
        j=i
        i=i+1
        if j%7==0:
            yield j

k=putNumbers(10)
print k gg
