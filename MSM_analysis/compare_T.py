import numpy as np
import argparse
from define_flow import define_flow_direction
from define_flow import read_cut
from define_flow import pos_pot

def calc_av(F, T, p):
	ms = np.shape(T)[0]

	return np.real( sum( sum( p[i] * T[i,j] *F[i,j] for j in range(ms)) for i in range(ms)  ))

def analyse_MSM(k):	
	Ev,Evec = np.linalg.eig(np.transpose(k))
	idx = Ev.argsort()[::-1] 
	Ev = Ev[idx]  #order eigenvectors by size
	Evec = np.transpose(Evec)
	Evec = Evec[idx,:] 
	return Ev, Evec



parser = argparse.ArgumentParser()
parser.add_argument('-o', help='ident')
parser.add_argument('-c', help='count')
parser.add_argument('-l', help='lagtime')
args = parser.parse_args()
start = int(args.o)
count = int(args.c)
lags = int(args.l)


ms = 30

T = np.zeros((count,ms,ms))

ident = str(start)
ms,cut = read_cut(ident)
r = define_flow_direction(ms,cut)

J = np.zeros(count)

for k in range(1,lags):
	for i in range(count):
		ident = str(start + i *10)
		T[i] = np.loadtxt("/data/isilon/bause/single_particle/T/T_"+str(k)+"_"+ident +".dat")

		Ev,Evec =  analyse_MSM(T[i])
		p = Evec[0] / sum(Evec[0]) 
		J[i] = calc_av(r,T[i],p);
		
	Jav = sum(J) / count 
	print(k, Jav)




