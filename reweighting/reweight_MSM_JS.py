import math
import argparse
import numpy as np
import os
from define_flow import define_flow_direction
from define_flow import read_cut
from define_flow import read_minpos
from define_flow import pos_pot
#from pyemma.msm import MSM

from define_cross import define_cross
from get_pos import get_pos
from fpt_ana import  first_passage#(start, end, T, length):
from fct_markov import mfpt_markov
from fct_markov import calc_Sprod
from fct_markov import analyse_MSM
from fct_markov import hist_to_cumhist
from fct_markov import calc_Caliber
from fct_markov import define_pathE
from fct_markov import calc_av
from fct_markov import moment_distribution

def reweight(T,r,sp, gamma,nu, ms):
	#W = [ [ np.power( T[i,j] , 2./3.) * np.power(T[j,i], 1./3.) * math.exp(2./3.*gamma*r[i,j] + 1./3. *gamma*r[j,i]) for j in range(ms)] for i in range(ms)]
	W = [ [  T[i,j] * math.exp(gamma*r[i,j] + nu* sp[i,j]) for j in range(ms)] for i in range(ms)]
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
	p = np.absolute(np.asarray(p)/sum(p))
	k = [[ Evecr[0,j] / (Evecr[0,i] * Ev[0]) *W[i,j] for j in range(ms)] for i in range(ms)]
	k= np.asarray(k)

	for i in range(ms):
		k[i] = k[i] / sum(k[i,:])

	return np.real(k),np.real(p)



parser = argparse.ArgumentParser()
parser.add_argument('-o', help='reweight to ident')
parser.add_argument('-i', help='reweight from ident')
parser.add_argument('-l', help='choose lagtime of T')
parser.add_argument('-rc', help='rancut for size of FPT well')
args = parser.parse_args();
ident = args.o
identin = args.i
rancut = int(args.rc)
lag = args.l


path = "/data/isilon/bause/single_particle/"+str(lag)+"/T/T_"+str(identin) +".dat"
T = np.loadtxt(path)
ms = np.shape(T)[0]
minT = np.nanmin(T)
for i in range(ms):
	for j in range(ms):
		if (T[i,j] > 0. ):
			T[i,j] = T[i,j]
		else:
			T[i,j] = minT

path0 = "/data/isilon/bause/"
ms,cut = read_cut(path0,ident)
ms, minpos = read_minpos(path0,ident)
r = define_flow_direction(ms,cut)


#T = np.loadtxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/T/Topt_"+str(identin)+".dat")
#T_comp = np.loadtxt("/data/isilon/bause/single_particle/"+str(lag)+"/T/T_"+str(ident)+".dat")

for line in open("/data/isilon/bause/single_particle/param_"+str(ident)+".dat","r"):
	cont = line.split()
	if(cont[0] == 'dT'):
		dT = float(cont[1])
	if(cont[0] == 'T0'):
		temperature = int(cont[1])
	if(cont[0] == 'extf'):
		force = int(cont[1])



potential = np.loadtxt("/data/isilon/bause/single_particle/potential_ms_"+str(ident)+".dat")

path = "/data/isilon/bause/single_particle/"+str(lag)+"/T/T_"+str(ident) +".dat"
T_comp = np.loadtxt(path)
minT = np.nanmin(T_comp)
for i in range(ms):
	for j in range(ms):
		if (T_comp[i,j] > 0. ):
			T_comp[i,j] = T_comp[i,j]
		else:
			T_comp[i,j] = minT


Ev,Evec =  analyse_MSM(T_comp)
p_comp = Evec[0] / sum(Evec[0]) 
temp = np.loadtxt("/data/isilon/bause/single_particle/potential_ms_"+str(ident)+".dat")
E = np.zeros((ms,ms))

for line in open("/data/isilon/bause/single_particle/param_"+str(ident)+".dat","r"):
	cont = line.split()
	if(cont[0] == 'dT'):
		dT = float(cont[1])
	if(cont[0] == 'extf'):
		force = int(cont[1])
	if(cont[0] == 'T'):
		temperature = float(cont[1])



potential = np.loadtxt("/data/isilon/bause/single_particle/potential_ms_"+str(ident)+".dat")

for i in range(ms):
	for j in range(ms):
		E[i,j] = (temp[i,1] + temp[j,1]) / 2.


lvl = len(minpos)
Ng = 100
gamma = np.zeros(Ng)
nu = np.zeros(Ng)
MFPT_ar = np.zeros((Ng, lvl *lvl))
ratio = np.zeros((Ng,3))

J_ar = np.zeros(Ng)
S_ar = np.zeros(Ng)

MFPT_comp = mfpt_markov(T_comp , minpos, rancut, ms,lvl).reshape(lvl*lvl)
ratio_comp =  sum( [MFPT_comp[1]/MFPT_comp[3] , MFPT_comp[6]/MFPT_comp[2] , MFPT_comp[5]/MFPT_comp[7]] ) / 3.
J_comp =  calc_av(r, T_comp, p_comp)


sp = np.zeros((ms,ms))
for k in range(ms):
	for l in range(ms):
		if (np.abs(k-l) < (ms/2)):
				sp[k,l] = ( -force / ms * (k-l)+ potential[k,1] - potential[l,1] ) / temperature
		elif (l>k):
				sp[k,l] = ( -force / ms * (ms-l+k)+ potential[k,1] - potential[l,1] ) / temperature
		else:
				sp[k,l] = ( force / ms * (ms-k+l)+ potential[k,1] - potential[l,1] ) / temperature


S_comp =  calc_av(E, T_comp, p_comp)
Sprod_comp = calc_Sprod(T_comp)
Eav_comp = calc_av(E,T_comp,p_comp)

gamma_start = -0.1
nu_start = -1
gamma_cur = gamma_start
nu_cur = nu_start
#ratio_cur = 10
J_cur = 0
err = 100


for i in range(1):
	for j in range(Ng):
		gamma[i] =   i*2*abs(gamma_start)/Ng + gamma_start  #random choice here
		gamma[i] = 0.
		nu[j] =   j*2*abs(nu_start)/Ng + nu_start 
		Trew, prew = reweight(T, r,sp, gamma[i],nu[j] , ms)
		J_ar[i] = calc_av(r,Trew,prew)
		S_ar[j] = calc_av(sp,Trew,prew)
		#Cerr = np.abs(J_ar[i] - J_comp) + np.abs(S_ar[j] - S_comp) 
		Cerr = np.abs(S_ar[j] - S_comp) 
		#print(gamma[i], nu[j] , S_ar[j])
		if ( Cerr < err):
			J_err = np.abs(J_ar[i] - J_comp) + np.abs(S_ar[j] - S_comp)
			gamma_cur = gamma[i]
			nu_cur = nu[j]
			T_cur = Trew[:,:]
			p_cur = prew[:]
			J_cur = J_ar[i]
			S_cur = S_ar[j]
			err = Cerr
	
#print( Trew[0,1], T_cur[0,1])
#print(T_cur)
#print(Trew)

Sprod_rew = calc_Sprod(T_cur)
MFPT_rew = mfpt_markov(T_cur , minpos, rancut, ms,lvl).reshape(lvl*lvl)
ratio_rew =  np.sum([MFPT_rew[1]/MFPT_rew[3] , MFPT_rew[6]/MFPT_rew[2] , MFPT_rew[5]/MFPT_rew[7]]) / 3.
Eav_rew = calc_av(E, T_cur, p_cur)


Ev,Evec =  analyse_MSM(T)
p_or = Evec[0] / sum(Evec[0]) 


P_cur = np.transpose([np.real(p_cur)])
P_comp = np.transpose([np.real(p_comp)])
P_or = np.transpose([np.real(p_or)])

out = np.concatenate((P_cur, P_comp, P_or),axis = 1)


np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/T/p_J_"+str(identin)+"_"+str(ident)+".dat",out)


#ratio_comp = sum(np.abs([MFPT_comp[1]/MFPT_comp[3] , MFPT_comp[6]/MFPT_comp[2] , MFPT_comp[5]/MFPT_comp[7] ] ))
#print(calc_Caliber(T_cur,p_cur,T,J_comp,r,ms,gamma_cur))


length = 400
MFPT_hist = np.zeros((lvl*lvl*3,length))
mom_cur = np.zeros((lvl*lvl,3))
mom_comp = np.zeros((lvl*lvl,3))
mom_start = np.zeros((lvl*lvl,3))
for i in range(lvl):
	for j in range(lvl):
		if (i != j):
			start = list(range(minpos[i]-rancut, minpos[i]+rancut+1))
			end = list(range(minpos[j]-rancut, minpos[j]+rancut+1))
			MFPT_hist[i*lvl+j] = first_passage(start, end, T_cur, length)
			MFPT_hist[i*lvl+j+   lvl*lvl] = first_passage(start,end, T_comp, length)
			MFPT_hist[i*lvl+j+ 2*lvl*lvl] = first_passage(start,end, T, length)
			for k in range(1,4):
				mom_cur[i*lvl+j,k-1] = moment_distribution(k,MFPT_hist[i*lvl+j])
				mom_comp[i*lvl+j,k-1] = moment_distribution(k,MFPT_hist[i*lvl+j+lvl*lvl])
				mom_start[i*lvl+j,k-1] = moment_distribution(k,MFPT_hist[i*lvl+j+2*lvl*lvl])

for k in range(1,4):
	np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/MFPT/moment_J_"+str(int(k))+"_"+str(identin)+"_"+str(ident)+".dat",mom_cur[:,k-1].reshape((lvl,lvl)))
	np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/MFPT/moment_J_"+str(int(k))+"_"+str(ident)+".dat",mom_comp[:,k-1].reshape((lvl,lvl)))


print(identin, ident,gamma_cur, nu_cur, J_comp, J_cur, S_comp, S_cur )

ratio = np.zeros((ms,ms))
act = np.zeros((ms,ms))
cutoff = 0.000001
for i in range(ms):
	for j in range(ms):
		if (T[i,j]> cutoff):
			ratio[i,j] = T_cur[i,j] / T_cur[j,i]
		act[i,j] = T_cur[i,j] * T_cur[j,i]

#print(identin, ident,gamma_cur, tot_error)
np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/MFPT/hist_J_"+str(identin)+"_"+str(ident)+".dat", np.transpose(MFPT_hist))
np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/T/T_J_"+str(identin)+"_"+str(ident)+".dat", T_cur)
np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/T/ratio_J_"+str(identin)+"_"+str(ident)+".dat",ratio)
np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/T/act_J_"+str(identin)+"_"+str(ident)+".dat", act)

	
gamma = np.transpose([gamma])
J_ar = np.transpose([J_ar])

out = np.concatenate((gamma, J_ar),axis = 1)

np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/J/J_"+str(identin)+".dat", out)




