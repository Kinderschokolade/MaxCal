import numpy as np

from msmtools.analysis import mfpt

def mfpt_markov(L, minpos, rancut, ms,lvl):
	mfpt_markov = np.zeros((lvl,lvl))
	for k in range(lvl):
		target = list(range(minpos[k] - rancut , minpos[k] + rancut))
		for t in range(lvl):
			origin = list(range(minpos[t] - rancut , minpos[t] + rancut))
			#print(origin, target)
			if(origin[-1] >= ms or target[-1] >= ms):
				print ("origin/target for mfpt out of range!") 
			else:
				mfpt_markov[t,k] = mfpt(L,target = target, origin = origin)
	return mfpt_markov


def calc_Sprod(T):
	ms = T.shape[1]
	S= 0.

	Ev,Evec = np.linalg.eig(np.transpose(T)) # trasnspose? 
	idx = Ev.argsort()[::-1] 
	Ev = Ev[idx]  #order eigenvectors by size
	Evec = np.transpose(Evec)
	Evec = Evec[idx,:] 
	p = np.real(Evec[0]  / sum(Evec[0])) #normalise

	for i in range(ms):
		for j in range(ms):
			if (T[i,j] > 0. and T[j,i] > 0. ):
				S+= p[i] *T[i,j] *np.log(T[i,j]/T[j,i])
	return S		


def analyse_MSM(k):	
	Ev,Evec = np.linalg.eig(np.transpose(k))
	idx = Ev.argsort()[::-1] 
	Ev = Ev[idx]  #order eigenvectors by size
	Evec = np.transpose(Evec)
	Evec = Evec[idx,:] 
	return Ev, Evec


def hist_to_cumhist(hist):
	shape = np.shape(hist)
	N = shape[0]
	length = shape[1]
	#print(length)
	cumhist = np.zeros((N,length))
	cumhist[:,0] = hist[:,0]
	for i in range(1,length):
		cumhist[:,i] = cumhist[:,i-1] + hist[:,i]
	return cumhist

def calc_Caliber(T,p,q,J,F,ms,gamma):
	S =0;
	for i in range(ms):
		for j in range(ms):
			if (not np.isclose(T[i,j],0.)):
				S+= -p[i]*T[i,j]* np.log(T[i,j]/q[i,j]) + gamma * p[i] *T[i,j] *F[i,j]
	
	S-= J

	return S

def define_pathE(ms,cut, ident): 
	r = np.zeros((ms,ms))
	lvl = np.size(cut)		
	E = np.loadtxt("/data/isilon/bause/single_particle/potential_ms_"+str(ident)+".dat")
	for i in range(ms):
		for j in range(ms):
			r[i,j] = (E[i] + E[j]) / 2.

	#np.savetxt("./Uij.dat",r)
	return r


def define_deltaE(ms,cut, ident, identin): 
	r = np.zeros((ms,ms))
	lvl = np.size(cut)
		
	temp = np.loadtxt("/data/isilon/bause/single_particle/potential_ms_"+str(ident)+".dat")
	tempin = np.loadtxt("/data/isilon/bause/single_particle/potential_ms_"+str(identin)+".dat")
	dE = tempin[:,1] - temp[:,1] 
	for i in range(ms):
		for j in range(ms):
			r[i,j] = (dE[i] + dE[j]) / 2.

	#np.savetxt("./Uij.dat",r)
	return r


def calc_av(F, T, p):
	ms = np.shape(T)[0]

	return np.real( sum( sum( p[i] * T[i,j] *F[i,j] for j in range(ms)) for i in range(ms)  ))

def moment_distribution(mom,dist):
	out =0
	for i in range(len(dist)):
		out += i**mom * dist[i]

	return out
	
