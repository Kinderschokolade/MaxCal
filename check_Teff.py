import numpy as np
import subprocess
import math
import argparse
from  define_flow import minpos
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


parser = argparse.ArgumentParser()
parser.add_argument('-o','--output', help='ident of output')
parser.add_argument('-T1', help='Temperature 1')
parser.add_argument('-T2', help='Temperature 2')
parser.add_argument('-T3', help='Temperature 3')
parser.add_argument('-U1', help='Potential 1')
parser.add_argument('-U2', help='Potential 2')
parser.add_argument('-U3', help='Potential 3')
parser.add_argument('-gamma', help='bath coupling constant')
parser.add_argument('-dT', help='time step size')
parser.add_argument('-fullT', help='simulation time')
parser.add_argument('-k', help='potential sharpnes')
args = parser.parse_args();
T1 = args.T1
T2 = args.T2
T3 = args.T3
U1 = args.U1
U2 = args.U2
U3 = args.U3
gamma = args.gamma
filename = args.output
fullT = args.fullT
dT = args.dT
potk = args.k
path = '/data/isilon/bause/single_particle/'

for i in range(-5,5):
	U1vary= str( int(U1) + i*0.5)
	filename = args.output+str(i+5)
	cmd = "~/code/single_particle/a.out -U1 "+U1vary+" -U2 "+U2+" -U3 "+U3+" -T1 "+T1+" -T2 "+T2+" -T3 "+T3+" -gamma "+gamma+" -o "+filename+" -fullT "+fullT+" -pot 2 -k "+potk +" -dT "+ dT
	sim = subprocess.Popen(cmd, shell=True)
	sim.wait()
	
cmd = 	"grep logl "+path+"res_"+args.output+"* > "+path +"log_"+args.output+".dat"
grepstep=subprocess.Popen(cmd,shell = True ,stdout=subprocess.PIPE)
grepstep.wait()
grep = "Tval"


dU = load_column(path+"log_"+args.output+".dat", 2)
fraclog = load_column(path+"log_"+args.output+".dat",1)
m1,b = np.polyfit(fraclog,dU,1)

dU = load_column(path+"log_"+args.output+".dat", 4)
fraclog = load_column(path+"log_"+args.output+".dat",3)
m2,b = np.polyfit(fraclog,dU,1)

Teff = str((-m1-m2) / 2.)
print(T1, m1, m2)


#alpha = []
#for i in range(-5,5):  
#	U1vary= str( int(U1) + i*0.5)
#	filename = str( int(args.output) + 1 )+str(i+5)
#	sim = subprocess.Popen(["~/code/single_particle/a.out -U1 "+U1vary+" -U2 "+U2+" -U3 "+U3+" -T1 "+Teff+" -T2 "+Teff+" -T3 "+Teff+" -gamma "+gamma+" -o "+filename+" -fullT "+fullT+" -pot 2 -k 100 -dT 0.01"], shell=True, stdout=subprocess.PIPE)
#	sim.wait()
#	
#	neqfile = args.output + str(i+5)
#	grep = "Tval"
#	Tneq = load_row(path+"res_"+neqfile+".dat", grep) 
#	Teq = load_row(path+"res_"+filename+".dat",grep)
#	alpha.append(Tneq-Teq)
#
#alpha= np.asarray(alpha)
#dU = np.asarray(dU)
#savedat = np.column_stack((dU,alpha))
#datafile = open('alpha.dat','wb')
#np.savetxt(datafile, savedat)
#
	



