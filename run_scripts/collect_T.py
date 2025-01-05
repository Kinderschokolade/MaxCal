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

numT = 15
numU = 11
starto1 = 20
starto2 = 0
startT= 15

path = "/data/isilon/bause/single_particle/res_"

tran = np.zeros((11,numT * numU))
for i in range(0,numT):
	T = startT - i -1
	o1 = str(starto1 + 2 * i)
	for j in range(0,numU):
		o2 = str( starto2 + j )
		grep = "Tval"
		if j<10:
			filename = path + o1 + "0" + o2 + ".dat"
		else:
			filename = path+ o1 + o2 + ".dat"
		out = load_row(filename,grep)
		dU =  j * 0.5
		tran[0][i*numU+j] = T
		tran[1][i*numU+j] = dU
		for k in range(0,9):
			tran[k+2][i*numU+j] = out[k] 
		
np.savetxt('data/fullT_02.dat', tran.transpose())	
