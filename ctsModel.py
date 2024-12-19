import numpy as np
import argparse
from msmbuilder.msm import ContinuousTimeMSM
from msmbuilder.msm import MarkovStateModel
from msmtools.analysis import mfpt
from define_flow import read_minpos
from define_flow import read_cut
from scipy.linalg import expm

ident = '0001'

ms,minpos = read_minpos(ident)
ms,cut = read_cut(ident)

F = open("/data/isilon/bause/single_particle/count_"+ ident +".dat","r") 
line = F.readline()
lagtime = line.split()
lagtime.pop(0)
lagtime  = [float(i) for i in lagtime]

cols = len(lagtime)

source = open("/data/pckr146/bause/single_particle/trajectory_"+ident+".dat", "r")
num = 1000000
traj = []

for i in range(num):
	stuff = source.readline()
	traj.append(int(float(stuff.split()[0]) * ms) )


lagtime  = 60
cutoff = 10#choose............. 

MSM = ContinuousTimeMSM(lag_time = lagtime, ergodic_cutoff = cutoff,n_timescales =8 )
MSM2 = MarkovStateModel(lag_time = lagtime , ergodic_cutoff = cutoff,n_timescales =8 )
#MSM.summarize()
#MSM.countmat_ = count
MSM.fit(traj)
MSM2.fit(traj)

rancut = 3 #random choice

N = len(minpos)
mfpt_ar = np.zeros((N,N))	
mfpt_ar2 = np.zeros((N,N))	

 
for k in range(N):
	target = list(range(minpos[k] - rancut , minpos[k] + rancut))
	for t in range(N):
		origin = list(range(minpos[t] - rancut , minpos[t] + rancut))
		if(origin[-1] >= ms or target[-1] >= ms):
			print ("origin/target for mfpt out of range!") 
		else:
			mfpt_ar[t][k] = mfpt(MSM.ratemat_,target = target, origin = origin)
			mfpt_ar2[t][k] = mfpt(MSM.transmat_,target = target, origin = origin) 


#print(mfpt_ar)
print(mfpt_ar2)

#mat = np.zeros((ms,ms))
#lag = np.zeros(10)
#for i in range(1,11):
#	mat = expm(MSM.ratemat_ * 1./(i))	
#	Ev,Evec = np.linalg.eig(mat)
#	idx = Ev.argsort()[::-1] 
#	Ev = Ev[idx]  #order eigenvectors by size3
#	Evec = np.transpose(Evec)
#	Evec = Evec[idx,:] 
#	print(1./i ,mfpt(mat,target = [4,5,6,7,8,9], origin = [14,15,16,17,18,19]) / (i),-(1./i)/np.log( Ev[1:9])) 
	
	
	
		
#print(MSM.summarize())
np.savetxt("test.dat" , MSM.transmat_)
np.savetxt("test2.dat" , MSM2.transmat_)
#print(MSM.ratemat_)
#print(MSM.n_states_)
print(MSM.timescales_ )
print(MSM2.timescales_ )
#print(mfpt_ar)
#print(mfpt_ar2)
	













