import math
import argparse
import numpy as np
import os
from define_flow import define_flow_direction
from define_flow import read_cut
from define_flow import read_minpos
from define_flow import pos_pot
from pyemma.msm import MSM
from msmtools.analysis import mfpt
from define_cross import define_cross
from get_pos import get_pos
from fpt_ana import  first_passage#(start, end, T, length):

def mfpt_markov(L,p, minpos, rancut, ms,lvl):
	mfpt_markov = np.zeros((lvl,lvl))
	for k in range(lvl):
		target = list(range(minpos[k] - rancut , minpos[k] + rancut))
		for t in range(lvl):
			origin = list(range(minpos[t] - rancut , minpos[t] + rancut))
			#print(origin, target)
			if(origin[-1] >= ms or target[-1] >= ms):
				print ("origin/target for mfpt out of range!") 
			else:
				mfpt_markov[t,k] = mfpt(L,target = target, origin = origin, mu = p)
	return mfpt_markov


def calc_Sprod(T):
	ms = T.shape[1]
	S= 0.

	Ev,Evec = np.linalg.eig(np.transpose(T)) # trasnspose? 
	idx = Ev.argsort()[::-1] 
	Ev = Ev[idx]  #order eigenvectors by size
	Evec = np.transpose(Evec)
	Evec = Evec[idx,:] 
	p = Evec[0]  / sum(Evec[0]) #normalise
	for i in range(ms):
		for j in range(ms):
			if (T[i,j] > 0. and T[j,i] > 0. ):
				S+= p[i] * T[i,j] * math.log(T[i,j] / T[j,i] )	
	return np.real(S)		

def calc_Caliber(T,p,q,J,F,ms,gamma):
	S =0;
	for i in range(ms):
		for j in range(ms):
			if (not np.isclose(T[i,j],0.)):
				S+= -p[i]*T[i,j]* np.log(T[i,j]/q[i,j]) 
	
	S+= np.abs(calc_av(F,T,p) - J)

	return S

def reweight(T,r, gamma, ms):
	W = [ [ T[i,j] * math.exp(gamma*r[i,j]) for j in range(ms)] for i in range(ms)]
	W = np.asarray(W)

	Ev,Evecl = np.linalg.eig(np.transpose(W))
	idx = Ev.argsort()[::-1] 
	Ev = Ev[idx]  #order eigenvectors by size
	Evecl = np.transpose(Evecl)
	Evecl = Evecl[idx,:] 
	Ev,Evecr = np.linalg.eig(W)
	idx = Ev.argsort()[::-1] 
	Ev = Ev[idx]  #order eigenvectors by size
	Evecr = np.transpose(Evecr)
	Evecr = Evecr[idx,:] 
	p  = Evecl[0,:] * Evecr[0,:]
	p = np.asarray(p)
	if (p[0] < 0.):
		p =  -1. *p
	k = [[ Evecr[0,j] / (Evecr[0,i] * Ev[0]) *W[i,j] for j in range(ms)] for i in range(ms)]
	k= np.asarray(k)


	return np.real(k),np.real(p)


def analyse_MSM(k):	
	Ev,Evec = np.linalg.eig(np.transpose(k))
	idx = Ev.argsort()[::-1] 
	Ev = Ev[idx]  #order eigenvectors by size
	Evec = np.transpose(Evec)
	Evec = Evec[idx,:] 
	return Ev, Evec


def calc_av(F, T, p):
	ms = np.shape(T)[0]
	return np.real( sum( sum( p[i] * T[i,j] *F[i,j] for j in range(ms)) for i in range(ms)  ))




#def flux(p, k):
#	ms = p.shape[0]
#	# this cannot work, because Tpl and Tmi is unknown after MSM reweighting!!
#	return (sum(sum( p[i] * k[i,j] for j in range(ms)) for i in range(ms)))

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


parser = argparse.ArgumentParser()
parser.add_argument('-o', help='reweight to ident')
parser.add_argument('-i', help='reweight from ident')
parser.add_argument('-rc', help='cut for minimum area definition')
parser.add_argument('-l', help='choose lagtime of T')
args = parser.parse_args();
ident = args.o
identin = args.i
rancut = int(args. rc) 
lag = args.l




path = "/data/isilon/bause/single_particle/"+str(lag)+"/T/Tpl_"+str(identin) +".dat"
Tpl = np.loadtxt(path)
path = "/data/isilon/bause/single_particle/"+str(lag)+"/T/Tmi_"+str(identin) +".dat"
Tmi = np.loadtxt(path)
ms = np.shape(Tpl)[0]
for i in range(ms):
	for j in range(ms):
		if ((Tmi[i,j] + Tpl[i,j]) > 0. ):
			Tmi[i,j] = Tmi[i,j] / sum(Tmi[i,:]+Tpl[i,:])
			Tpl[i,j] = Tpl[i,j] / sum(Tpl[i,:]+Tmi[i,:])
		else:
			Tmi[i,j] = 0.
			Tpl[i,j] = 0.

T = Tmi+Tpl

path = "/data/isilon/bause/single_particle/"+str(lag)+"/T/Tpl_"+ str(ident) +".dat"
Tpl = np.loadtxt(path)
path = "/data/isilon/bause/single_particle/"+str(lag)+"/T/Tmi_"+ str(ident) +".dat"
Tmi = np.loadtxt(path)
ms = np.shape(Tpl)[0]
for i in range(ms):
	for j in range(ms):
		if ((Tmi[i,j] + Tpl[i,j]) > 0. ):
			Tmi[i,j] = Tmi[i,j] / sum(Tmi[i,:]+Tpl[i,:])
			Tpl[i,j] = Tpl[i,j] / sum(Tpl[i,:]+Tmi[i,:])
		else:
			Tmi[i,j] = 0.
			Tpl[i,j] = 0.

T_comp = Tmi+Tpl
Ev, Evec = analyse_MSM(T_comp)

p_comp = np.real(Evec[0] / sum(Evec[0,:]))

path="/data/isilon/bause/"
ms,cut = read_cut(path,ident)
ms, minpos = read_minpos(path,ident)
r = define_flow_direction(ms,cut)

np.savetxt("Fij.dat", r)
lvl = len(minpos)
Ng = 100
gamma = np.zeros(Ng)
MFPT_ar = np.zeros((Ng, lvl *lvl))
ratio = np.zeros((Ng,3))
J = np.zeros(Ng)

MFPT_comp = mfpt_markov(T_comp ,p_comp, minpos, rancut, ms,lvl).reshape(lvl*lvl)
ratio_comp =  [MFPT_comp[1]/MFPT_comp[3] , MFPT_comp[6]/MFPT_comp[2] , MFPT_comp[5]/MFPT_comp[7]]
J_comp = calc_av(r,T_comp,p_comp)

Sprod_comp = calc_Sprod(T_comp)

gamma_start = -0.1
gamma_cur = gamma_start
ratio_cur = 10

for i in range(Ng):
	gamma[i] =   i*2*abs(gamma_start)/Ng + gamma_start  #random choice here
	Trew, prew = reweight(T, r, gamma[i] , ms)
	Ev, Evec = analyse_MSM(Trew)
	MFPT_ar[i] = mfpt_markov(Trew,prew, minpos, rancut, ms,lvl).reshape(lvl*lvl)
	ratio[i] = [MFPT_ar[i][1]/MFPT_ar[i][3] , MFPT_ar[i][6]/MFPT_ar[i][2] , MFPT_ar[i][5]/MFPT_ar[i][7] ] 
	J[i] = calc_av(r,Trew,prew)
#	print(gamma[i] , flux(r,Trew,prew), MFPT_ar[i][1], MFPT_ar[i][3] , sum(ratio[i]))
	# for three-well case
	print(gamma[i],calc_Caliber(Trew,prew,T,J_comp,r,ms,gamma[i]))
	if ( sum(np.abs(ratio[i] - ratio_comp)) < ratio_cur):
		ratio_cur = sum(np.abs(ratio[i] - ratio_comp))
		ratcur = sum(np.abs(ratio[i]))
		gamma_cur = gamma[i]
		T_cur = Trew
		p_cur = prew
		MFPT_rew = MFPT_ar[i][1]
		J_cur = calc_av(r,Trew,prew)
		#print( Trew[0,1], T_cur[0,1])
#print(T_cur)
#print(Trew)

Sprod_rew = calc_Sprod(T_cur)

ratio_comp = sum(np.abs([MFPT_comp[1]/MFPT_comp[3] , MFPT_comp[6]/MFPT_comp[2] , MFPT_comp[5]/MFPT_comp[7] ] ))

print(calc_Caliber(T_cur,p_cur,T,J_comp,r,ms,gamma_cur))
print(identin, ident,gamma_cur,  MFPT_comp[1], MFPT_rew, J_comp, J_cur, Sprod_comp, Sprod_rew, ratio_comp, ratcur )

length = 1000
MFPT_hist = np.zeros((lvl*lvl*3,length))
for i in range(lvl):
	for j in range(lvl):
		if (i != j):
			start = list(range(minpos[i]-rancut, minpos[i]+rancut+1))
			end = list(range(minpos[j]-rancut, minpos[j]+rancut+1))
			MFPT_hist[i*lvl+j] = first_passage(start, end, T_cur, length)
			MFPT_hist[i*lvl+j+   lvl*lvl] = first_passage(start,end, T_comp, length)
			MFPT_hist[i*lvl+j+ 2*lvl*lvl] = first_passage(start,end, T, length)


tot_error= sum(sum(abs(MFPT_hist[i,j] - MFPT_hist[i+lvl*lvl,j]) for i in range(lvl*lvl)) for j in range(length))

MFPT_cumhist = hist_to_cumhist(MFPT_hist)
#print(identin, ident,gamma_cur, tot_error)
np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/MFPT/MFPThist_"+str(identin)+"_"+str(ident)+".dat", np.transpose(MFPT_hist))
np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/MFPT/MFPTcumhist_"+str(identin)+"_"+str(ident)+".dat", np.transpose(MFPT_cumhist))
np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/T/MFPTT_"+str(identin)+"_"+str(ident)+".dat", T_cur)

	
gamma = np.transpose([gamma])
J = np.transpose([J])

out = np.concatenate((gamma, ratio,J),axis = 1)

np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/ratio_"+str(identin)+".dat", out)



