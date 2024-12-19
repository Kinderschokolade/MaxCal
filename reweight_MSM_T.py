import math
import argparse
import numpy as np
import os
import getpass
from define_flow import define_flow_direction
from define_flow import read_cut
from define_flow import read_minpos
from define_flow import pos_pot
from pyemma.msm import MSM
from msmtools.analysis import mfpt
from define_cross import define_cross
from get_pos import get_pos
from fpt_ana import  first_passage#(start, end, T, length):

def define_pathE(ms,cut, ident, identin): 
	r = np.zeros((ms,ms))
	lvl = np.size(cut)
		
	temp = np.loadtxt("/data/isilon/bause/single_particle/potential_ms_"+str(ident)+".dat")
	#tempin = np.loadtxt("/data/isilon/bause/single_particle/potential_ms_"+str(identin)+".dat")
	#dE = tempin[:,1] - temp[:,1] 
	dE = temp[:,1]
	for i in range(ms):
		for j in range(ms):
			r[i,j] = (dE[i] + dE[j]) / 2.

	np.savetxt("./Uij.dat",r)
	return r

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
	p = Evec[0]  / sum(Evec[0]) #normalise
	for i in range(ms):
		for j in range(ms):
			if (T[i,j] > 0. and T[j,i] > 0. ):
				S+= p[i] * T[i,j] * math.log(T[i,j] / T[j,i] )	
	return np.real(S)		

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

	k = [[ Evecr[0,j] / (Evecr[0,i] * Ev[0]) *W[i,j] for j in range(ms)] for i in range(ms)]
	k= np.asarray(k)

	for i in range(ms):
		k[i] = k[i] / sum(k[i,:])
	p = p / sum(p)
	return np.real(k),np.real(p)


def analyse_MSM(k):	
	Ev,Evec = np.linalg.eig(np.transpose(k))
	idx = Ev.argsort()[::-1] 
	Ev = Ev[idx]  #order eigenvectors by size
	Evec = np.transpose(Evec)
	Evec = Evec[idx,:] 
	return Ev, Evec


def flux(p, k):
	ms = p.shape[0]
	# this cannot work, because Tpl and Tmi is unknown after MSM reweighting!!
	return (sum(sum( p[i] * k[i,j] for j in range(ms)) for i in range(ms)))

def E(p,k,E):
	ms = p.shape[0]
	return (np.real(sum( sum( p[i] * k[i,j] * E[i,j] for j in range(ms)) for i in range(ms))))


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

user = getpass.getuser()

if (user=="mbause"):
	path0 = "/u/mbause/data/"
	path1 = "/u/mbause/data/"
else:
	path0 = "/data/isilon/bause/"
	path1 = "/data/pckr194/bause/"






path = path0+"single_particle/T/Tpl_"+str(lag)+"_"+str(identin) +".dat"
Tpl = np.loadtxt(path)
path = path0+"single_particle/T/Tmi_"+str(lag)+"_"+str(identin) +".dat"
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

path = path0+"single_particle/T/Tpl_"+str(lag)+"_"+ str(ident) +".dat"
Tpl = np.loadtxt(path)
path = path0+"single_particle/T/Tmi_"+str(lag) +"_"+ str(ident) +".dat"
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

Ev,Evec =  analyse_MSM(T_comp)
p_comp = np.real(Evec[0] / sum(Evec[0]) )

ms,cut = read_cut(path0,ident)
ms, minpos = read_minpos(path0,ident)
#r = define_flow_direction(ms,cut)

r = define_pathE(ms,cut,ident,identin)

lvl = len(minpos)
Ng = 200
gamma = np.zeros(Ng)
MFPT_ar = np.zeros((Ng, lvl *lvl))
ratio = np.zeros((Ng,3))

MFPT_comp = mfpt_markov(T_comp , minpos, rancut, ms,lvl).reshape(lvl*lvl)
#ratio_comp =  [MFPT_comp[1]/MFPT_comp[3] , MFPT_comp[6]/MFPT_comp[2] , MFPT_comp[5]/MFPT_comp[7]]
Sprod_comp = calc_Sprod(T_comp)

temp = np.loadtxt(path0+"single_particle/potential_ms_"+str(ident)+".dat")
E_vec = np.zeros((ms,ms))
for i in range(ms):
	for j in range(ms):
		E_vec[i,j] =  (temp[i,1] + temp[j,1])/ 2. 
	
avE_comp = E(p_comp,T_comp,E_vec)
gamma_start = -1
gamma_cur = gamma_start
MFPT_cur = 10000
delta_E = 1000
for i in range(Ng):
	gamma[i] =   i*2*abs(gamma_start)/Ng + gamma_start  #random choice here
	Trew, prew = reweight(T, r, gamma[i] , ms)
	#Ev, Evec = analyse_MSM(Trew)
	MFPT_ar[i] = mfpt_markov(Trew, minpos, rancut, ms,lvl).reshape(lvl*lvl)
#	ratio[i] = [MFPT_ar[i][1]/MFPT_ar[i][3] , MFPT_ar[i][6]/MFPT_ar[i][2] , MFPT_ar[i][5]/MFPT_ar[i][7] ] 
	# for three-well case
	#np.savetxt("temp.dat", Trew)
	#np.savetxt("diff.dat", T_comp - Trew)
	avE = E(prew,Trew,E_vec)
	Sprod_rew = calc_Sprod(Trew)
#	print(gamma[i] ,avE_comp, avE, Sprod_comp, Sprod_rew, MFPT_ar[i,1] , MFPT_comp[1],  MFPT_ar[i,2] , MFPT_comp[2],  MFPT_ar[i,7] , MFPT_comp[7], MFPT_ar[i,6], MFPT_comp[6]) 
#	print(gamma[i] , MFPT_ar[i,1] , MFPT_comp[1],  MFPT_ar[i,2] , MFPT_comp[2],  MFPT_ar[i,7] , MFPT_comp[7], MFPT_ar[i,6], MFPT_comp[6]) 

	if ( sum(np.abs(MFPT_ar[i] - MFPT_comp)) < MFPT_cur):
#	if (np.abs(avE- avE_comp) < delta_E):
		delta_E = np.abs(avE-avE_comp)
		MFPT_cur = sum( np.abs( MFPT_ar[i] - MFPT_comp) )
		gamma_cur = gamma[i]
		p_cur = prew[:]
		T_cur = Trew[:,:]
		MFPT_rew = MFPT_ar[i]
#print(T_cur)
#print(Trew)

Sprod_rew = calc_Sprod(T_cur)
avE_cur = E(p_cur,T_cur,E_vec)

ratio_comp = sum(np.abs([MFPT_comp[1]/MFPT_comp[3] , MFPT_comp[6]/MFPT_comp[2] , MFPT_comp[5]/MFPT_comp[7] ] ))

print(identin, ident,gamma_cur, avE_comp, avE_cur,  MFPT_comp[1], MFPT_rew[1], MFPT_comp[2], MFPT_rew[2],MFPT_comp[5], MFPT_rew[5], Sprod_comp, Sprod_rew, ratio_comp, sum(np.abs([MFPT_ar[i,1]/MFPT_ar[i,3] , MFPT_ar[i,6]/MFPT_ar[i,2] , MFPT_ar[i,5]/MFPT_ar[i,7] ] )) )
#print(identin, ident, gamma_cur)
#print(avE_comp, avE_cur)
#print(MFPT_comp)
#print(MFPT_rew )
detbal = np.zeros((ms,ms))
detbal_comp = np.zeros((ms,ms))

#for i in range(ms):
#	for j in range(ms):
#		detbal[i,j] = p_cur[i] * T_cur[i,j] - p_cur[j] * T_cur[j,i]
#		detbal_comp[i,j] = p_comp[i] * T_comp[i,j] - p_comp[j] * T_comp[j,i]

#np.savetxt("detbal.dat", detbal)
#np.savetxt("detbal_comp.dat", detbal_comp)

np.savetxt("p1.dat", p_comp)
np.savetxt("p2.dat", p_cur)
#print(np.real(p_comp))
#print(p_cur)


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
np.savetxt("/data/isilon/bause/single_particle/reweight/MFPThist_"+str(identin)+"_"+str(ident)+".dat", np.transpose(MFPT_hist))
np.savetxt("/data/isilon/bause/single_particle/reweight/MFPTcumhist_"+str(identin)+"_"+str(ident)+".dat", np.transpose(MFPT_cumhist))
np.savetxt("/data/isilon/bause/single_particle/reweight/T_"+str(identin)+"_"+str(ident)+".dat", T_cur)

	
gamma = np.transpose([gamma])

out = np.concatenate((gamma, ratio),axis = 1)

np.savetxt("/data/isilon/bause/single_particle/reweight/ratio_"+str(identin)+".dat", out)



