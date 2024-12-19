import numpy as np
import subprocess
import math
import argparse
from time import sleep
import multiprocessing as mp
from functools import partial


parser = argparse.ArgumentParser()
parser.add_argument('-o','--output', help='ident of output')
parser.add_argument('-T1', help='Temperature 1')
parser.add_argument('-T2', help='Temperature 2')
parser.add_argument('-T3', help='Temperature 3')
parser.add_argument('-U1', help='Potential 1')
parser.add_argument('-U2', help='Potential 2')
parser.add_argument('-U3', help='Potential 3')
parser.add_argument('-gamma', help='bath coupling constant')
parser.add_argument('-dT', help='simulation time step')
parser.add_argument('-fullT', help='simulation time')
parser.add_argument('-k', help='potential sharpnes')
parser.add_argument('-f', help='external force')
parser.add_argument('-ms', help='microstates')
args = parser.parse_args();
T1 = args.T1
T2 = args.T2
T3 = args.T3
U1 = args.U1
U2 = args.U2
U3 = args.U3
f_ext = args.f
gamma = args.gamma
ms = args.ms
filename = args.output
fullT = args.fullT
dT = args.dT
potk = args.k
path = '/data/isilon/bause/single_particle/'

offset =0
steps = 10
stepsize = 0.5
endU = float(U1)+ steps*stepsize +0.5
startU = float(U1) - steps*stepsize +0.5


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
	#out.pop(0)
	fil =list(filter(lambda a: a != grep.encode() , out))
	fil = np.asarray(fil, dtype=np.float32)
	return(fil)

def run_Sim(U1,U2,U3,T1,T2,T3,gamma,filename,fullT,potk,dT,offset,ms):
	cmd = "~/code/single_particle/a.out -U1 "+str(U1)+" -U2 "+U2+" -U3 "+U3+" -T1 "+T1+" -T2 "+T2+" -T3 "+T3+" -gamma "+gamma+" -o "+filename+" -fullT "+fullT+" -pot 2 -k "+potk+" -dT "+dT + " -lvl 3 -ms " + str(ms)
	sim = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
	sim.wait()

# wrapper for run_Sim
def run_sim_wrap(i):
	U1vary = startU + i * stepsize
	if i<10: 
		filename = str( int(args.output) + offset )+"0"+str(i)
	else:
		filename = str( int(args.output) + offset )+str(i)

	return run_Sim(U1vary,U2,U3,T1,T2,T3,gamma,filename,fullT,potk,dT,offset,ms)


offset = 0
if __name__=='__main__':
	pool = mp.Pool()
	run_sim_u1_range = range(steps)
	pool.map_async(run_sim_wrap, run_sim_u1_range).get( 9999999 ) 	


mp.Barrier(pool)

	
cmd = 	"grep logl "+path+"res_"+args.output+"* > "+path +"log_"+args.output+".dat"
grepstep=subprocess.Popen(cmd,shell = True ,stdout=subprocess.PIPE)
grepstep.wait()
grep = "Tval"

# only valid for this special case with ms = 3
dU = load_column(path+"log_"+args.output+".dat", 2)
fraclog = load_column(path+"log_"+args.output+".dat",1)
m1,b = np.polyfit(fraclog,dU,1)

dU = load_column(path+"log_"+args.output+".dat", 4)
fraclog = load_column(path+"log_"+args.output+".dat",3)
m2,b = np.polyfit(fraclog,dU,1)

Teff = str((-m1-m2) / 2.)
print("noneq-sim",-m1, -m2)

alpha = []
newline = ["\n"]

offset=1
T1 = Teff
T2 = Teff
T3 = Teff
if __name__=='__main__':
	pool = mp.Pool()
	run_sim_u2_range = range(steps)
	pool.map_async(run_sim_wrap, run_sim_u2_range).get( 9999999 ) 	
	
mp.Barrier(pool)
	
eqpos= str(int(args.output) +1 )
cmd = 	"grep logl "+path+"res_"+eqpos+"* > "+path +"log_"+eqpos+".dat"
grepstep=subprocess.Popen(cmd,shell = True ,stdout=subprocess.PIPE)
grepstep.wait()
grep = "Tval"

dU = load_column(path+"log_"+eqpos+".dat", 2)
fraclog = load_column(path+"log_"+eqpos+".dat",1)
m1,b = np.polyfit(fraclog,dU,1)

dU = load_column(path+"log_"+eqpos+".dat", 4)
fraclog = load_column(path+"log_"+eqpos+".dat",3)
m2,b = np.polyfit(fraclog,dU,1)

print("eq-sim",-m1, -m2)


for i in range(steps):
	if i<10: 
		neqfile = str( int(args.output) )+"0"+str(i)
		filename =  str( int(args.output) +1 )+"0"+str(i)
	else:
		neqfile = str( int(args.output) )+str(i)
		filename =  str( int(args.output) +1 )+str(i)

	grep = "Tval"
	Tneq = load_row(path+"res_"+neqfile+".dat", grep) 
	Teq = load_row(path+"res_"+filename+".dat",grep)
	grep= "Terr"
	Tneq_err = load_row(path+"res_"+neqfile+".dat", grep)
	Teq_err = load_row(path+"res_"+filename+".dat", grep)
	out = np.concatenate((Tneq-Teq,Tneq_err - Teq_err),axis=0)
	alpha.append(out)
	
alpha=np.column_stack((dU,alpha))


open('data/alpha_'+args.output+'.dat', 'wb').close() #clear file
datafile = open('data/alpha_'+args.output+'.dat','ab')

np.savetxt(datafile,alpha)

eqout = str( int(args.output) +1 )	
cmd = 	"grep logl "+path+"res_"+eqout+"* > "+path +"log_"+eqout+".dat"
grepstep=subprocess.Popen(cmd,shell = True ,stdout=subprocess.PIPE)
grepstep.wait()
	
dU = load_column(path+"log_"+eqout+".dat", 4)
fraclog = load_column(path+"log_"+eqout+".dat",3)
m2,b = np.polyfit(fraclog,dU,1)

dU = load_column(path+"log_"+eqout+".dat", 2)
fraclog = load_column(path+"log_"+eqout+".dat",1)
m1,b = np.polyfit(fraclog,dU,1)


Teff = str((-m1-m2) / 2.)
print("eq-sim",-m1, -m2)


