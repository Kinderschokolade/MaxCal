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
parser.add_argument('-l', help='write T of specific line')
parser.add_argument('-rc', help='rancut', default=1)
args = parser.parse_args();
ident = args.o
lag = int(args.l)
rancut = int(args.rc)

ms, minpos = read_minpos(ident)
ms, cut = read_cut(ident)

L = np.zeros((ms,ms))
lvl = len(minpos)
mfpt_ar = np.zeros((lvl,lvl))

blocks = 100

path = "/data/pckr146/bause/single_particle/trajectory_" +str(ident) +".dat"
series = open(path, "r")
tra_len = int(os.popen("wc -l "+path).readline().split()[0])-blocks
block_len = int(tra_len / (blocks*lag))

print(tra_len)

stuff =series.readline()
split = stuff.split()
Gamma_old = float(split[0])
pos_old = get_pos(Gamma_old,ms)

Count = np.zeros((ms,ms))
q = np.zeros(ms)	
length = 1000
fpt_ana = np.zeros((blocks,length))
mfpt_ana = np.zeros(blocks)	


# successively increase number of blocks to check convergence
for b in range(blocks):
	#read information
	for i in range(block_len):
		for j in range(lag-1):
			series.readline()

		stuff = series.readline()
		split = stuff.split()
		Gamma = float(split[0])
		v = float(split[1])
		pos = get_pos(Gamma,ms)
		Count[pos_old, pos] +=1
		q[pos] +=1
		pos_old = pos

	for k in range(ms):
		for l in range(ms):
			if (q[k] > 0.):
				L[k,l] = Count[k,l] / q[k]	
			else:
				L[k,l] = 0.
				L[k,k] = 1. 
	#find p stationary
	Ev,Evec = np.linalg.eig(np.transpose(L))
	idx = Ev.argsort()[::-1] 
	Ev = Ev[idx]  #order eigenvectors by size3
	Evec = np.transpose(Evec)
	Evec = Evec[idx,:] 
	p_stat = Evec[0]
		
	# calc fpt (only 1->2)
	start = list(range(minpos[0]-rancut, minpos[0]+rancut+1))
	end = list(range(minpos[1]-rancut, minpos[1]+rancut+1))
	fpt_ana[b,:] = first_passage(start,end,L,length+1)
	mfpt_ana[b] = sum( i * fpt_ana[b,i] for i in range(length))
	

np.savetxt("/data/isilon/bause/single_particle/FPT/hist_compare_trajectory_"+str(lag)+"_"+ident+".dat", np.transpose(fpt_ana))
np.savetxt("/data/isilon/bause/single_particle/FPT/mfpt_compare_trajectory_"+str(lag)+"_"+ident+".dat", mfpt_ana)

##vary micostates

series.seek(0)
mstates = [12,18,24,30,36,42,48,51,60,90,120,150,180,210,240,270,300,360,420,480,540,600,720,840]

numms = len(mstates)

Count = []
q = []
L = []
p_stat = []
minpos = []

rmin = [0.25,0.25+ 1./3., 0.25 + 2./3.]

pos_old = np.zeros(numms)
pos = np.zeros(numms)

fpt_ana = np.zeros((numms,length))
mfpt_ana  = np.zeros((2,numms))


for i in range(numms):
	Count.append(np.zeros((mstates[i],mstates[i])))
	q.append(np.zeros(mstates[i]))
	L.append(np.zeros((mstates[i],mstates[i])))
	p_stat.append(np.zeros(mstates[i]))
	minpos.append( [get_pos(rmin[0],mstates[i]), get_pos(rmin[1],mstates[i]), get_pos(rmin[2],mstates[i]) ] )


tra_len = int(tra_len /lag)
for i in range(tra_len):
		for j in range(lag-1):
			series.readline()

		stuff = series.readline()
		split = stuff.split()
		Gamma = float(split[0])
		for j in range(numms):
			pos[j] = get_pos(Gamma,mstates[j])
			Count[j][pos_old[j], pos[j]] +=1
			q[j][pos[j]] +=1
			pos_old[:] = pos



for i in range(numms):
	for k in range(mstates[i]):
		for l in range(mstates[i]):
			if (q[i][k] > 0.):
				L[i][k,l] = Count[i][k,l] / q[i][k]
			else:
				L[i][k,l] = 0.
				L[i][k,k] = 1. 
	#find p stationary
	Ev,Evec = np.linalg.eig(np.transpose(L[i]))
	idx = Ev.argsort()[::-1] 
	Ev = Ev[idx]  #order eigenvectors by size3
	Evec = np.transpose(Evec)
	Evec = Evec[idx,:] 
	p_stat[i] = Evec[0]
		
	# calc fpt (only 1->2)
	rancut =  int( mstates[i] / 30)	
	print(mstates[i] , (minpos[i][0] -rancut)/mstates[i] ,minpos[i][0]/mstates[i], (minpos[i][0] +rancut +1)/mstates[i] )

	start = list(range(minpos[i][0]-rancut, minpos[i][0]+rancut+1))
	end = list(range(minpos[i][1]-rancut, minpos[i][1]+rancut+1))

	fpt_ana[i,:] = first_passage(start,end,L[i],length+1)
	mfpt_ana[1,i] = sum( k * fpt_ana[i,k] for k in range(length))
	
mfpt_ana[0] = mstates


np.savetxt("/data/isilon/bause/single_particle/FPT/hist_compare_ms_"+str(lag)+"_"+ident+".dat", np.transpose(fpt_ana))
np.savetxt("/data/isilon/bause/single_particle/FPT/mfpt_compare_ms_"+str(lag)+"_"+ident+".dat", np.transpose(mfpt_ana))





















	
