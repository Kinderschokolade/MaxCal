import numpy as np
import subprocess
import math
import argparse

def load_column(filename, column):
	p = []
	i=0
	with open(filename) as inf:
		for line in inf:
			parts = line.split() # split line into parts
			p.append(float(parts[column]))
	p=np.array(p)
	return(p)


def load_row (filename, grep):
	thing =subprocess.Popen(["grep "+grep+" "+filename],shell=True, stdout=subprocess.PIPE).communicate()[0]
	out = thing.split()
	out.pop(0)
	out = np.asarray(out, dtype=np.float32)
	return(out)

num = 10
numU = 10
starto = 40
startT= 5


alpha = np.zeros((20,numU))
outfile = open('data/fullalpha_01.dat','ab')

for i in range(0,num):
	o = str(starto + 2 *i)
	filename = '/home/theorie/bause/code/single_particle/data/alpha_'+o+'.dat'
	T = startT + i
	dU = load_column(filename,0)
	for j in range(0,numU): 
		alpha[0][j] = T
		alpha[1][j] = dU[j]

	for k in range(0,18):
		stuff = load_column(filename, k+1)
		for j in range(0,numU):
			alpha[k+2][j] = stuff[j] 
	np.savetxt(outfile,alpha.transpose())

	np.savetxt(outfile,[''],fmt="%s")
		
