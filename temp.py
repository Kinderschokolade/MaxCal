import subprocess
import math
import argparse
import numpy as np
from time import sleep


def load_row (filename, grep):
	thing =subprocess.Popen(["grep "+grep+" "+filename],shell=True, stdout=subprocess.PIPE).communicate()[0]
	out = thing.split()
	out.pop(0)
	out = np.asarray(out, dtype=np.float32)
	return(out)



parser = argparse.ArgumentParser()
parser.add_argument('-o','--output', help='ident of output')
args = parser.parse_args();

alpha = []
dU = np.zeros(20);

path = '/data/isilon/bause/single_particle/'
for i in range(-10,10):  
	dU[i+10] = 5 + i *0.5 
	if i<0: 
		filename = str( int(args.output) + 1 )+"0"+str(i+10)
	else:
		filename = str( int(args.output) + 1 )+str(i+10)

	if (i<0):
		neqfile = args.output +"0"+ str(i+10)
	else:
		neqfile = args.output + str(i+10)

	grep = "Tval"
	Tneq = load_row(path+"res_"+neqfile+".dat", grep) 
	Teq = load_row(path+"res_"+filename+".dat",grep)
	print(Tneq)
	grep= "Terr"
	Tneq_err = load_row(path+"res_"+neqfile+".dat", grep)
	Teq_err = load_row(path+"res_"+filename+".dat", grep)
	out = np.concatenate((Tneq-Teq,Tneq_err - Teq_err),axis=0)
	alpha.append(out)

alpha= np.asarray(alpha)

savedat = np.column_stack((dU,alpha))
datafile = open('data/alpha_'+args.output+'.dat','wb')
np.savetxt(datafile, savedat)

