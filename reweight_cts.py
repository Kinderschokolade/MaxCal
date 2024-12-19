import math
import argparse
import numpy as np
import os
from define_flow import define_flow_direction
from define_flow import read_cut
from define_flow import read_minpos
from define_flow import pos_pot
#from pyemma.msm import MSM
from msmtools.analysis import mfpt
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
from scipy.optimize import fsolve

def equ(c,ms,W):
	tup = [ (sum([W[j,k] *c[k]*c[j] for k in range(ms)] )-1.0) for j in range(ms)]
	tup = tuple(tup)
	return tup


def reweight(gamma,T,r,E, ms,J_comp, Sprod, Sprod_in):

	W = [ [  T[i,j] * math.exp(gamma/2.*(r[i,j]+r[j,i]) + (Sprod[i,j]-Sprod_in[i,j])/2. ) for j in range(ms)] for i in range(ms)]
	W = np.asarray(np.real(W))


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

	c = np.real(Evecr[0,:])

	copt = fsolve(equ,c,args=(ms,W,),xtol=1e-8) # should converge for any starting point

	k = [[copt[j] *copt[i] *W[i,j] for j in range(ms)] for i in range(ms)]
	k= np.asarray(k)

	for i in range(ms): # minimal correction from numerics, eta = 1.
		k[i,j] = k[i,j] / sum(k[i,:]) 
	
	Ev,Evec =  analyse_MSM(k)
	p = np.real(Evec[0] / sum(Evec[0]) )
	J = calc_av(r,k,p)
	J_err =  np.abs(J - J_comp)

	return J_err, np.real(k),np.real(p)



parser = argparse.ArgumentParser()
parser.add_argument('-i', help='reweight from ident')
parser.add_argument('-l', help='choose lagtime of T')
parser.add_argument('-rc', help='rancut for size of FPT well')
args = parser.parse_args();
identin = args.i
rancut = int(args.rc)
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

path0 = "/data/isilon/bause/"
ms,cut = read_cut(path0,identin)
ms, minpos = read_minpos(path0,identin)
r = define_flow_direction(ms,cut)


#T = np.loadtxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/T/Topt_"+str(identin)+".dat")
#T_comp = np.loadtxt("/data/isilon/bause/single_particle/"+str(lag)+"/T/T_"+str(ident)+".dat")

for line in open("/data/isilon/bause/single_particle/param_"+str(identin)+".dat","r"):
	cont = line.split()
	if(cont[0] == 'dT'):
		dT = float(cont[1])
	if(cont[0] == 'T0'):
		temperature = int(cont[1])
	if(cont[0] == 'extf'):
		force = int(cont[1])



potential = np.loadtxt("/data/isilon/bause/single_particle/potential_ms_"+str(identin)+".dat")
path = "/data/isilon/bause/single_particle/"+str(lag)+"/T/Tpl_"+ str(identin) +".dat"
Tpl = np.loadtxt(path)
path = "/data/isilon/bause/single_particle/"+str(lag)+"/T/Tmi_"+ str(identin) +".dat"
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

#############3
#path = "/data/isilon/bause/single_particle/"+str(lag)+"/reweight/T/Topt_3020.dat"
#T_comp = np.loadtxt(path)
#############3

Ev,Evec =  analyse_MSM(T_comp)
p_comp = Evec[0] / sum(Evec[0]) 
temp = np.loadtxt("/data/isilon/bause/single_particle/potential_ms_"+str(identin)+".dat")
E = np.zeros((ms,ms))

for i in range(ms):
	for j in range(ms):
		E[i,j] = (temp[i,1] + temp[j,1]) / 2.

potential = temp [:,1]

lvl = len(minpos)
Ng = 200
gamma = np.zeros(Ng)

f_start = -10
f_cur = f_start
delta_f = 2*(-f_start)/Ng

J_comp = calc_av(r,T_comp,p_comp)

length = 400
mom1 = np.zeros((lvl*lvl+1,Ng))
mom2 = np.zeros((lvl*lvl+1,Ng))
mom3 = np.zeros((lvl*lvl+1,Ng))
Sprod = np.zeros((2,Ng))

MFPT_hist = np.zeros((lvl*lvl, length))

p_all = np.zeros((Ng,ms+1))

for line in open("/data/isilon/bause/single_particle/param_"+str(identin)+".dat","r"):
	cont = line.split()
	if(cont[0] == 'dT'):
		dT = float(cont[1])
	if(cont[0] == 'T0'):
		temperature = int(cont[1])
	if(cont[0] == 'extf'):
		force_in = int(cont[1])


Sprod_out = np.zeros((ms,ms))
Sprod_in = np.zeros((ms,ms))
for k in range(ms):
	for l in range(ms):
		if (np.abs(k-l) <15):
				Sprod_in[k,l] = ( -force_in / ms * (k-l)+ potential[k] - potential[l] ) / temperature
		elif (l>k):
				Sprod_in[k,l] = ( -force_in / ms * (ms-l+k)+ potential[k] - potential[l] ) / temperature
		else:
				Sprod_in[k,l] = ( force_in / ms * (ms-k+l)+ potential[k] - potential[l] ) / temperature


for i in range(Ng):
	gamma = 0. #does not matter for now
	f_cur = f_start + i * delta_f
	for k in range(ms):
		for l in range(ms):
			if (np.abs(k-l) <15):
					Sprod_out[k,l] = ( -f_cur / ms * (k-l)+ potential[k] - potential[l] ) / temperature
			elif (l>k):
					Sprod_out[k,l] = ( -f_cur / ms * (ms-l+k)+ potential[k] - potential[l] ) / temperature
			else:
					Sprod_out[k,l] = ( f_cur / ms * (ms-k+l)+ potential[k] - potential[l] ) / temperature


	J_err,Trew,prew=reweight(gamma,T_comp,r,E, ms,J_comp, Sprod_out, Sprod_in)
	J = calc_av(r,Trew,prew)
	Sp = calc_Sprod(Trew)

	p_all[i,0] = J	
	p_all[i,1:] = prew

	mom1[0,i] = J
	mom2[0,i] = J
	mom3[0,i] = J
	Sprod[0,i] = J
	Sprod[1,i] = Sp

	for l in range(lvl):
		for m in range(lvl):
			if (i != j):
				start = list(range(minpos[l]-rancut, minpos[l]+rancut+1))
				end = list(range(minpos[m]-rancut, minpos[m]+rancut+1))
				MFPT_hist[l*lvl+m] = first_passage(start, end, Trew, length)
				mom1[l*lvl+m+1,i] = moment_distribution(1,MFPT_hist[l*lvl+m])
				mom2[l*lvl+m+1,i] = moment_distribution(2,MFPT_hist[l*lvl+m])
				mom3[l*lvl+m+1,i] = moment_distribution(3,MFPT_hist[l*lvl+m])


np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/T/p_cts_"+str(identin) +".dat",p_all)

np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/MFPT/mom1_cts_"+str(identin) +".dat",np.transpose(mom1))
np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/MFPT/mom2_cts_"+str(identin) +".dat",np.transpose(mom2))
np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/MFPT/mom3_cts_"+str(identin) +".dat",np.transpose(mom3))
np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/MFPT/Sprod_cts_"+str(identin) +".dat",np.transpose(Sprod))


