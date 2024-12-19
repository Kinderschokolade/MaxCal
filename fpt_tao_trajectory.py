import argparse
import numpy as np
import math
import os
from define_flow import read_minpos
from define_flow import read_cut
from fpt_ana import first_passage

def get_pos(x,ms):
	pos = math.floor(x*ms)
	if(pos == ms):
		 pos -= 1
	return pos



parser = argparse.ArgumentParser()
parser.add_argument('-o', help='ident')
parser.add_argument('-ms', help='write T of specific line')
parser.add_argument('-rc', help='rancut', default=1)
parser.add_argument('-l', help='maximal lagtime', default=1)
args = parser.parse_args();
ident = args.o
maxlag = int(args.l)
rancut = int(args.rc)
ms = int(args.ms)

path = "/data/pckr146/bause/single_particle/trajectory_" +str(ident) +".dat"
series = open(path, "r")
tra_len = int(os.popen("wc -l "+path).readline().split()[0])-maxlag

print(tra_len)

pos_old = np.zeros(maxlag)
stuff =series.readline()
split = stuff.split()
Gamma_old = float(split[0])
pos_old[:] = get_pos(Gamma_old,ms)

Count = np.zeros((ms,ms))
q = np.zeros(ms)	
length = 1000
fpt_ana = np.zeros((maxlag,length))
mfpt_ana = np.zeros(maxlag)	

##vary lagtime

Count = np.zeros((maxlag,ms,ms))
q = np.zeros((maxlag,ms))
L = np.zeros((maxlag,ms,ms))
p_stat = np.zeros((maxlag,ms))
rmin = [0.25,0.25+ 1./3., 0.25 + 2./3.]
minpos = [get_pos(rmin[0],ms), get_pos(rmin[1],ms), get_pos(rmin[2],ms) ] 

fpt_ana = np.zeros((maxlag,length))
fpt_ana_back = np.zeros((maxlag,length))
mfpt_ana  = np.zeros((2,maxlag))
mfpt_ana_back  = np.zeros((2,maxlag))

for i in range(tra_len):
	stuff = series.readline()
	split = stuff.split()
	Gamma = float(split[0])
	pos = get_pos(Gamma,ms)
	for j in range(1,maxlag):
		if (i % j == 0):
			Count[j,pos_old[j], pos] +=1
			q[j][pos] +=1
			pos_old[j] = pos



for i in range(maxlag):
	for k in range(ms):
		for l in range(ms):
			if (q[i][k] > 0.):
				L[i,k,l] = Count[i,k,l] / q[i,k]
			else:
				L[i,k,l] = 0.
				L[i,k,k] = 1. 
	#find p stationary
	Ev,Evec = np.linalg.eig(np.transpose(L[i]))
	idx = Ev.argsort()[::-1] 
	Ev = Ev[idx]  #order eigenvectors by size3
	Evec = np.transpose(Evec)
	Evec = Evec[idx,:] 
	p_stat[i] = Evec[0]
		
	# calc fpt (only 1->2)

	start = list(range(minpos[0]-rancut, minpos[0]+rancut+1))
	end = list(range(minpos[1]-rancut, minpos[1]+rancut+1))

	fpt_ana[i,:] = first_passage(start,end,L[i],length)
	mfpt_ana[1,i] = sum( k * fpt_ana[i,k] for k in range(length))

	start = list(range(minpos[1]-rancut, minpos[1]+rancut+1))
	end = list(range(minpos[0]-rancut, minpos[0]+rancut+1))
	fpt_ana_back[i,:] = first_passage(start,end,L[i],length+1)
	mfpt_ana_back[1,i] = sum( k * fpt_ana_back[i,k] for k in range(length))


mfpt_ana[0] = list(range(maxlag))
mfpt_ana_back[0] = list(range(maxlag))


np.savetxt("/data/isilon/bause/single_particle/FPT/hist_compare_tao_"+str(ms)+"_"+ident+".dat", np.transpose(fpt_ana))
np.savetxt("/data/isilon/bause/single_particle/FPT/mfpt_compare_tao_"+str(ms)+"_"+ident+".dat", np.transpose(mfpt_ana))

np.savetxt("/data/isilon/bause/single_particle/FPT/hist_compareb_tao_"+str(ms)+"_"+ident+".dat", np.transpose(fpt_ana_back))
np.savetxt("/data/isilon/bause/single_particle/FPT/mfpt_compareb_tao_"+str(ms)+"_"+ident+".dat", np.transpose(mfpt_ana_back))






















	
