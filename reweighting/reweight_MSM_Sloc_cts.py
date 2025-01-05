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

from scipy.linalg import eig
from scipy.optimize import fmin
from scipy.optimize import fsolve
from scipy.optimize import least_squares

import random

def equ(c,ms,W):
	tup = [ (sum([W[j,k] * c[k]*c[j] for k in range(ms)] )-1.0) for j in range(ms)]
	tup = tuple(tup)
	return tup

def reweight_Inv(Invariant,Sprod,ms):

#	W = [ [  T[i,j] * math.exp(gamma/2.*(r[i,j]+r[j,i]) + (Sprod[i,j]-Sprod_in[i,j])/2. ) for j in range(ms)] for i in range(ms)]
	W = [ [  Invariant[i,j]/1000. * math.exp( (Sprod[i,j])/2. ) for j in range(ms)] for i in range(ms)]

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

	copt = fsolve(equ,c,args=(ms,W,),xtol=1e-8) # should converge for any starting point
#	copt = least_squares(equ, c,args=(ms,W,), bounds=([0.]*ms, 1.))
	k = [[copt[j] *copt[i] *W[i,j] for j in range(ms)] for i in range(ms)]
	k= np.asarray(k)
	print(copt)
	for i in range(ms):
		norm = sum(k[i,:])
		for j in range(ms): # minimal correction from numerics, eta = 1.
			k[i,j] = k[i,j] / norm
		print(norm)
	
	Ev,Evec =  analyse_MSM(k)
	p = np.real(Evec[0] / sum(Evec[0]) )

	return np.real(k),np.real(p)

def reweight(T,Sprod,ms):

#	W = [ [  T[i,j] * math.exp(gamma/2.*(r[i,j]+r[j,i]) + (Sprod[i,j]-Sprod_in[i,j])/2. ) for j in range(ms)] for i in range(ms)]
	W = [ [  np.sqrt(T[i,j] * T[j,i]) * math.exp( (Sprod[i,j])/2. ) for j in range(ms)] for i in range(ms)]

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

#	c = np.abs(np.real(Evecr[0,:]))
	c = np.ones(ms)
#	copt = fsolve(equ,c,args=(ms,W,),xtol=1e-4) # should converge for any starting point
	copt = least_squares(equ, c,args=(ms,W,), bounds=([0.]*ms, 1.))
	copt = copt.x
	k = [[copt[j] *copt[i] *W[i,j] for j in range(ms)] for i in range(ms)]
	k= np.asarray(k)
	for i in range(ms):
		norm = sum(k[i,:])
		for j in range(ms): # minimal correction from numerics, eta = 1.
			k[i,j] = k[i,j] / norm

	Ev,Evec =  analyse_MSM(k)
	p = np.real(Evec[0] / sum(Evec[0]) )

	return np.real(k),np.real(p)






parser = argparse.ArgumentParser()
parser.add_argument('-o', help='reweight to ident')
parser.add_argument('-l', help='choose lagtime of T')
parser.add_argument('-rc', help='rancut for size of FPT well')
parser.add_argument('-fm', help='maximal force')
parser.add_argument('-s', help='numer of steps -fm')

args = parser.parse_args();
ident = args.o
rancut = int(args.rc)
lag = args.l
force_max = float(args.fm)
steps = int(args.s)


path = "/data/isilon/bause/single_particle/"+str(lag)+"/T/T_"+str(ident) +".dat"
T = np.loadtxt(path)
ms = np.shape(T)[0]
for i in range(ms):
	for j in range(ms):
		if (T[i,j]  > 0. ):
			T[i,j] = T[i,j] / sum(T[i,:])
		else:
			T[i,j] = 0.


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
	if(cont[0] == 'T'):
		temperature_in = float(cont[1])
	if(cont[0] == 'extf'):
		force_in = int(cont[1])
	if(cont[0] == 'gamma'):
		gamma = int(cont[1])


potential = np.loadtxt("/data/isilon/bause/single_particle/potential_ms_"+str(ident)+".dat")


temp = np.loadtxt("/data/isilon/bause/single_particle/potential_ms_"+str(ident)+".dat")
E = np.zeros((ms,ms))
 

for i in range(ms):
	for j in range(ms):
		E[i,j] = (temp[i,1] + temp[j,1]) / 2.

minpos= [3,20,36,52]

lvl = len(minpos)

#ratio_cur = 10
J_cur = 0
J_err = 100

Jdiff = 1.
Jdiff_min = 0.00001

sp = np.zeros((ms,ms))
sp_in = np.zeros((ms,ms))

#for k in range(ms):
#	for l in range(ms):
#		if (np.abs(k-l) < (ms/2)):
#				sp_in[k,l] = ( -force_in / ms * (k-l)+ potential[k,1] - potential[l,1] ) / temperature_in
#		elif (l>k):
#				sp_in[k,l] = ( -force_in / ms * (ms-l+k)+ potential[k,1] - potential[l,1] ) / temperature_in
#		else:
#				sp_in[k,l] = (force_in / ms * (ms-k+l)+ potential[k,1] - potential[l,1] ) / temperature_in


################33

for k in range(ms):
	for l in range(ms):
		if (np.abs(k-l) < (ms/2)):
			dist = (k-l)/ms
		elif (l>k):
			dist = (ms-l+k)/ms
		else:
			dist = (ms -k +l)/ms
		
		dist = dist*dist

#		T[k,l] = np.sqrt(gamma/(temperature_in)) * np.exp(-dist*gamma/(4.*temperature_in*(float(lag))*dT))
#
#for i in range(ms):
#	norm = sum(T[i,:])
#	for j in range(ms):
#			T[i,j] = T[i,j] / norm
###########reweighting from free diffusion does not work well... What am I forgetting?



#Invariant = 2.*np.loadtxt("T.dat")
#path = "/data/isilon/bause/single_particle/"+str(lag)+"/reweight/T/Inv_opt_"+str(ident) +".dat"
path = "/data/isilon/bause/single_particle/"+str(lag)+"/T/T_"+str(ident) +".dat"
#Invariant = np.loadtxt(path)
T_ref = np.loadtxt(path)
#Invariant = np.transpose(Invariant)
for i in range(ms):
	for j in range(ms):
		if (T_ref[i,j] > 0. ):
			T_ref[i,j] = T_ref[i,j]
		else:
			T_ref[i,j] = 0.



###################

#print( "real",p_comp)
#gamma_cur = fmin(reweight_min,x0=0.0,ftol = Jdiff_min , disp=False, args=(T ,r,E , ms , J_comp, sp, sp_in))
#gamma_cur = gamma_cur[0]


#for i in range(3):
#	gamma_cur = i * -0.01
#	J_err,T_cur,p_cur = reweight(gamma_cur,T ,r , ms , J_comp, sp)
#	print(gamma_cur,J_err)

gamma_cur =0 

J_cts = np.zeros((2,steps))
Sprod_cts = np.zeros((2,steps))
MFPT_cts = np.zeros((1+lvl*lvl,steps))
MFPT2_cts = np.zeros((1+lvl*lvl,steps))
MFPT3_cts = np.zeros((1+lvl*lvl,steps))

length = 10000
J_comp = 0.
MFPT_hist = np.zeros((lvl*lvl,length))


p_hist = np.zeros((lvl+1,steps))	

sp_loc = np.zeros(ms)
for i in range(steps):

	f_loc = np.zeros(ms)

	f_loc[1]=  1.90367e-05  
	f_loc[2]=  0.000535806  
	f_loc[3]=  0.00819126  
	f_loc[4]=  0.0680176  
	f_loc[5]=  0.306775  
	f_loc[6]=  0.75153  
	f_loc[7]=  1  
	f_loc[8]=  0.722739 
	f_loc[9]=  0.283721 
	f_loc[10]=  0.0604961 
	f_loc[11]=  0.00700636 
	f_loc[12]=  0.000440742 
	f_loc[13]=  1.50593e-05 


	f_loc[:] = f_loc[:] * force_max/ms *i /steps

	#f_loc[9:17] = +force_max /ms * i/steps

	for k in range(ms-1):
		sp_loc[k] = (potential[k,1] - potential[k+1,1] + f_loc[k] ) /temperature_in

	sp_loc[ms-1] = (potential[ms-1,1] - potential[0,1] + f_loc[ms-1] ) /temperature_in


	for k in range(ms):
		for l in range(ms):
			if (np.abs(k-l) < (ms/2)):
				if (k<l):
					sp[k,l] = sum(sp_loc[k:l])
				else:
					sp[k,l] = -sum(sp_loc[l:k])
				#sp_in[k,l] = ( -force_in / ms * (k-l)*mul+ potential_in[k,1] - potential_in[l,1] ) / temperature_in
			elif (l>k):
				sp[k,l] = -sum(sp_loc[0:k]) - sum(sp_loc[l:ms])
				#sp[k,l] = ( -force / ms * (ms-l+k)*mul+ potential[k,1] - potential[l,1] ) / temperature
				#sp_in[k,l] = ( -force_in / ms * (ms-l+k)*mul+ potential_in[k,1] - potential_in[l,1] ) / temperature_in
			else:
				sp[k,l] = -sum(sp_loc[0:k]) -sum(sp_loc[l:ms])+ 1*(sum(f_loc[0:k]) + sum(f_loc[l:ms]))
				#sp[k,l] = ( force / ms * (ms-k+l)*mul+ potential[k,1] - potential[l,1] ) / temperature
				#sp_in[k,l] = (force_in / ms * (ms-k+l)*mul+ potential_in[k,1] - potential_in[l,1] ) / temperature_in
#	J_err, T_cur, p_cur = reweight(gamma_cur,T,r,E,ms,J_comp,sp,sp_in)
	T_cur, p_cur = reweight(T_ref,sp,ms)

#	J_cts[1,i] = calc_av(r,T_cur,p_cur)
  #	Sprod_cts[1,i] = calc_Sprod(T_cur)
	#ratio_rew =  np.sum([MFPT_rew[1]/MFPT_rew[3] , MFPT_rew[6]/MFPT_rew[2] , MFPT_rew[5]/MFPT_rew[7]]) / 3.
	#Eav_rew = calc_av(E, T_cur, p_cur)
	#p_cur # track lvl of minima?
	force = force_max * i/steps
	J_cts[0,i] = force;
	MFPT_cts[0,i] = force;
	MFPT2_cts[0,i] = force;
	MFPT3_cts[0,i] = force;
	Sprod_cts[0,i] = force;
	if (np.any(p_cur > 1.)):
		p_hist[:,i] = np.nan
	else:
		p_hist[0,i] = force
		for k in range(lvl):
			p_hist[k+1,i] = sum(p_cur[minpos[k]-rancut:minpos[k]+rancut+1])
			
	

#	print(p_hist[:,i])

	for k in range(lvl):
		for l in range(lvl):
			if (i != j):
				start = list(range(minpos[k]-rancut, minpos[k]+rancut+1))
				end = list(range(minpos[l]-rancut, minpos[l]+rancut+1))
				MFPT_hist[k*lvl+l] = first_passage(start, end, T_cur, length)
					
				MFPT_cts[k*lvl+l+1,i] = moment_distribution(1, MFPT_hist[k*lvl+l])
				MFPT2_cts[k*lvl+l+1,i] = moment_distribution(2,MFPT_hist[k*lvl+l])
				MFPT3_cts[k*lvl+l+1,i] = moment_distribution(3,MFPT_hist[k*lvl+l])

				MFPT2_cts[k*lvl+l+1,i] = np.sqrt( MFPT2_cts[k*lvl+l+1,i] - MFPT_cts[k*lvl+l+1,i]* MFPT_cts[k*lvl+l+1,i] ) #standard deviation
				MFPT3_cts[k*lvl+l+1,i] = (MFPT3_cts[k*lvl+l+1,i] - 3.*MFPT_cts[k*lvl+l+1,i] * MFPT2_cts[k*lvl+l+1,i] * MFPT2_cts[k*lvl+l+1,i] -  MFPT_cts[k*lvl+l+1,i]* MFPT_cts[k*lvl+l+1,i] * MFPT_cts[k*lvl+l+1,i] )/ (  MFPT2_cts[k*lvl+l+1,i] * MFPT2_cts[k*lvl+l+1,i] * MFPT2_cts[k*lvl+l+1,i])


np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/ps_cts_"+str(ident)+".dat",np.transpose(p_hist))

np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/mom1_cts_"+str(ident)+".dat",np.transpose(MFPT_cts))
np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/mom2_cts_"+str(ident)+".dat",np.transpose(MFPT2_cts))
np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/mom3_cts_"+str(ident)+".dat",np.transpose(MFPT3_cts))
#np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/J_cts_"+str(ident)+".dat",np.transpose(J_cts))
#np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/Sprod_cts_"+str(ident)+".dat",np.transpose(Sprod_cts))

