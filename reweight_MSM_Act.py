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

from scipy.misc import logsumexp
from scipy import optimize

import random

def equ(c,ms,W):
	tup = [ (sum([W[j,k] *c[k] * c[j] for k in range(ms)] )-1.0) for j in range(ms)]
	tup = tuple(tup)
	return tup


def reweight(T, ms, act_out):
	W = np.zeros((ms,ms))
	for i in range(ms):
		for j in range(ms):
			if ( T[j,i] > 0. ):	
				W[i,j] = math.sqrt(T[i,j] / T[j,i] * act_out[i,j] )
	
	np.savetxt("W.dat", W)

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
	copt = np.abs(np.real(Evecl[0,:]))
	#c = np.zeros(ms);
	k = [[copt[j] / (copt[i] * Ev[0]) *W[i,j] for j in range(ms)] for i in range(ms)]
	k= np.asarray(k)
	for i in range(ms): # minimal correction from numerics, eta = 1.
		norm = sum(k[i,:])
		k[i,:] = k[i,:] / norm

	Ev,Evec =  analyse_MSM(k)
	p = np.real(Evec[0] / sum(Evec[0]) )
	J = calc_av(r,k,p)
	J_err =  np.abs(J - J_comp)
	
	return  np.real(k),np.real(p),copt


def equ_v(v,ms,W,d_A):
	tup = [ (sum([W[i,j] * (d_A[i,k] + v[j] ) * np.exp((v[i]-v[j])/2.) for j in range(ms)] )- v[i]) for i in range(ms)]
	tup = tuple(tup)
	return tup



def reweight2(T, ms, act_out):
	W = np.zeros((ms,ms))
	d_A = np.zeros((ms,ms))
	for i in range(ms):
		for j in range(ms):
			if ( T[j,i] > 0. ):	
				W[i,j] = math.sqrt(T[i,j] / T[j,i] * act_out[i,j] )
				d_A[i,j] = np.log( T[i,j] * T[j,i] ) - act_out[i,j] 
	
	v = np.ones(ms)	
	v = np.optimize.fsolve(equ_v,v,args=(ms,W,d_A,),xtol=1e-12)

	k = [[  for j in range(ms)] for i in range(ms)]
	k= np.asarray(k)
	for i in range(ms): # minimal correction from numerics, eta = 1.
		norm = sum(k[i,:])
		k[i,:] = k[i,:] / norm

	Ev,Evec =  analyse_MSM(k)
	p = np.real(Evec[0] / sum(Evec[0]) )
	J = calc_av(r,k,p)
	J_err =  np.abs(J - J_comp)
	
	return  np.real(k),np.real(p),copt



ms_max = 120
p_glob = np.zeros((ms_max,ms_max)) 
alpha_glob = np.zeros(ms_max) 
term = np.zeros((ms_max,ms_max)) 
def minimiser(x,ms,W,w,Sprod,Sprod_in):
	global p_glob
	global alpha_glob
	for i in range(ms):
		for j in range(ms):
			p_glob[i,j] =W[i,j]*np.exp(w[j,i] *(x[ms+i]+x[j] )  +w[i,j] * (x[ms+j] + x[i]) )
		alpha_glob[i] = sum(  p_glob[i,:]*w[i,:] *(Sprod[i,:] -Sprod_in[i,:] -x[ms+i] +x[ms:] -x[:ms] + x[i]) )

	func = sum( (sum(p_glob[i,:] ) -1. )**2. + (x[ms+i] +x[i] +1. + alpha_glob[i])**2.   for i in range(ms)) #+ 0.0001 *np.linalg.norm(x)
	#print(func)
	return func

def callback(x, f, acc):	
	func1 =  1*sum( abs(sum(p_glob[i,:] ) - 1. ) for i in range(ms)) #+ 0.0001 *np.linalg.norm(x)
	func2 =  1*sum( abs(x[i] - x[2*ms] + alpha_glob[i]-1.) for i in range(ms)) #+ 0.0001 *np.linalg.norm(x)
	print(func1, func2)

def minimiser_trafo(x,ms,W,w,Sprod,Sprod_in):
	global p_glob
	global alpha_glob
	for i in range(ms):
		alpha_glob[i] = 0.
		for j in range(ms):
			p_glob[i,j] =W[i,j]*np.exp(w[i,j] *(-x[ms+i] + x[ms+j] )  + 1./2. * (x[i] + x[j] + x[ms+i] - x[ms+j]) )
			alpha_glob[i] += p_glob[i,j] *w[i,j] *(Sprod[i,j] - Sprod_in[i,j] - x[ms+i] + x[ms+j] )

	func1 =  sum( np.abs(sum(p_glob[i,:]) -1. ) for i in range(ms)) #+ 0.0001 *np.linalg.norm(x)
	func2 =  sum( np.abs(x[i] - x[2*ms] + alpha_glob[i]-1.) for i in range(ms)) #+ 0.0001 *np.linalg.norm(x)
	func = 2.*func1+func2
	return func

def reweight_full(T,ms,Sprod):
	Omega,k,p,copt = reweight(T, ms, Sprod)
	w = np.zeros((ms,ms))
	cutoff = 0.0001
	for i in range(ms):
		for j in range(ms):
			if (k[j,i] > cutoff):
				w[i,j] = 1./ (1+ p[i] / p[j] * np.exp(Sprod[i,j]))
			else:
				w[i,j] = 1./2.

	np.savetxt("w.dat", w)
			

	W = np.zeros((ms,ms))

	Sprod_in = np.zeros((ms,ms))
	for i in range(ms):
		for j in range(ms):
#			W[i,j] = T[i,j]**w[j,i] * T[j,i]**w[i,j] *np.exp(w[i,j]* Sprod[i,j] -1.)	
			if(T[i,j] > cutoff and T[j,i] > cutoff):
				Sprod_in[i,j] = np.log( T[i,j] / T[j,i] )
			W[i,j] = T[i,j] *np.exp( w[i,j] *(Sprod[i,j] - Sprod_in[i,j] ) - 1.)	

	x = np.zeros(2*ms+1)
	x[:ms] = copt[:]

	for ii in range(3):
		minimizer_kwargs = { "args": (ms,W,w,Sprod,Sprod_in),"method":"L-BFGS-B",}# "options":{"gtol": 1e-4,"ftol":1e-6}, }
	#	out = optimize.basinhopping(minimiser_trafo, x,minimizer_kwargs=minimizer_kwargs, callback = callback, disp = 'True', niter = 10, interval = 100 , stepsize=0.5 )
		out = optimize.minimize(minimiser_trafo,x,args=(ms,W,w,Sprod,Sprod_in), method="L-BFGS-B", options={"gtol":1e-02,"ftol":1e-6})
		x = out.x
		v = x[:ms]
		u = x[ms:2*ms]
		zeta =x[2*ms]
		k = [ [ W[i,j]*np.exp(w[i,j] *(-u[i]+u[j] )  +1./2.* (v[i] + v[j] - u[j] + u[i] )  ) for j in range(ms) ] for i in range(ms)]	
		k = np.asarray(k)
		knorm = np.zeros(ms)
		
		print("norm",sum(np.abs([sum( k[i,:]) - 1. for i in range(ms)]))  ) # should be 1.
		np.savetxt("norm.dat", [sum( k[i,:]) - 1. for i in range(ms)] ) # should be 1.
		print("glob", sum( np.abs([ -zeta+ +v[i] +sum( w[i,j] * k[i,j]*(Sprod[i,j] - Sprod_in[i,j]  -u[i] + u[j] )  for j in range(ms)) - 1. for i in range(ms) ])) ) # should be 1
		np.savetxt("glob.dat", [ -zeta+ +v[i] +sum( w[i,j] * k[i,j]*(Sprod[i,j] - Sprod_in[i,j]  -u[i] + u[j] )  for j in range(ms)) - 1. for i in range(ms) ]) # should be 1
	
	
		for i in range(ms):
			knorm[i] = sum(k[i,:])
			k[i,:] /= knorm[i]
	
		Invariant = np.zeros((ms,ms))
		norm = 0.
		m = 1./2.*(u+v)
		n = 1./2.*(v-u)
		for i in range(ms):
			for j in range(ms):
				Invariant[i,j] = T[i,j] * T[j,i] * np.exp(-1 + 1./2.*(v[i]+v[j] + u[i] - u[j] ) + w[i,j] * (Sprod[i,j] - Sprod_in[i,j] +u[j] - u[i] ) )
				#norm += Invariant[i,j]
	
	#	Invariant = Invariant / norm
	
	
		Ev,Evec =  analyse_MSM(k)
		pold = p
		p = np.real(Evec[0] / sum(Evec[0]) )
		print("diff",ii,sum(np.abs( pold- p)) )
		for i in range(ms):
			for j in range(ms):
				if (k[j,i] > cutoff):
					w[i,j] = 1./ (1+ p[i] / p[j] * np.exp(Sprod[i,j]))
				else:
					w[i,j] = 1./2.	
		for i in range(ms):
			for j in range(ms):
				W[i,j] = T[i,j] *np.exp( w[i,j] *(Sprod[i,j] - Sprod_in[i,j] ) - 1.)	


	return Invariant,np.real(k),np.real(p), x


parser = argparse.ArgumentParser()
parser.add_argument('-o', help='reweight to ident')
parser.add_argument('-i', help='reweight from ident')
parser.add_argument('-l', help='choose lagtime of T')
parser.add_argument('-rc', help='rancut for size of FPT well')
parser.add_argument('-a', help='rescaling for \delta x',default=1.0)
parser.add_argument('-d', help='exit condition for iteration',default=1e-05)
args = parser.parse_args();
ident = args.o
identin = args.i
rancut = int(args.rc)
lag = args.l
mul = float(args.a)
delta_min = float(args.d)


path = "/data/isilon/bause/single_particle/"+str(lag)+"/T/T_"+str(identin) +".dat"
T = np.loadtxt(path)

path0 = "/data/isilon/bause/"
ms,cut = read_cut(path0,ident)
ms, minpos = read_minpos(path0,ident)

for i in range(ms):	
	for j in range(ms):
		if (T[i,j] > 0.):	
			T[i,j] = T[i,j] 
		else:
		#	T[i,j] = minT
			T[i,j] = 0. # doesnt make a difference! minT might help to find data wher pij = 0


#------------------------
#minpos = np.asarray([4,19,35,52])
#---------------------


r = define_flow_direction(ms,cut)

#T = np.loadtxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/T/Topt_"+str(identin)+".dat")
#T_comp = np.loadtxt("/data/isilon/bause/single_particle/"+str(lag)+"/T/T_"+str(ident)+".dat")

for line in open("/data/isilon/bause/single_particle/param_"+str(ident)+".dat","r"):
	cont = line.split()
	if(cont[0] == 'dT'):
		dT = float(cont[1])
	if(cont[0] == 'T0'):
		temperature = float(cont[1])
	if(cont[0] == 'extf'):
		force = int(cont[1])
	if(cont[0] == 'T'):
		temperature = float(cont[1])

for line in open("/data/isilon/bause/single_particle/param_"+str(identin)+".dat","r"):
	cont = line.split()
	if(cont[0] == 'dT'):
		dT_in = float(cont[1])
	if(cont[0] == 'T0'):
		temperature_in = float(cont[1])
	if(cont[0] == 'extf'):
		force_in = int(cont[1])
	if(cont[0] == 'T'):
		temperature_in = float(cont[1])




potential = np.loadtxt("/data/isilon/bause/single_particle/potential_ms_"+str(ident)+".dat")
potential_in = np.loadtxt("/data/isilon/bause/single_particle/potential_ms_"+str(identin)+".dat")
path = "/data/isilon/bause/single_particle/"+str(lag)+"/T/T_"+ str(ident) +".dat"
T_comp = np.loadtxt(path)
ms = np.shape(T_comp)[0]
for i in range(ms):
	for j in range(ms):
		if (T_comp[i,j] > 0. ):
			T_comp[i,j] = T_comp[i,j]
		else:
			T_comp[i,j] = 0.

Ev,Evec =  analyse_MSM(T_comp)
p_comp = np.real(Evec[0] / sum(Evec[0]) )

Ev,Evec =  analyse_MSM(T)
p_in = np.real(Evec[0] / sum(Evec[0]) )

temp = np.loadtxt("/data/isilon/bause/single_particle/potential_ms_"+str(ident)+".dat")
E = np.zeros((ms,ms))
 

for i in range(ms):
	for j in range(ms):
		E[i,j] = np.exp(  (temp[i,1] + temp[j,1]) / (2. * temperature) )


lvl = len(minpos)
Ng = 800
gamma = np.zeros(Ng)
MFPT_ar = np.zeros((Ng, lvl *lvl))
ratio = np.zeros((Ng,3))

J_ar = np.zeros(Ng)

MFPT_comp = mfpt_markov(T_comp , minpos, rancut, ms,lvl).reshape(lvl*lvl)
#ratio_comp =  sum( [MFPT_comp[1]/MFPT_comp[3] , MFPT_comp[6]/MFPT_comp[2] , MFPT_comp[5]/MFPT_comp[7]] ) / 3.
J_comp =  calc_av(r, T_comp, p_comp)


Sprod_comp = calc_Sprod(T_comp)
Eav_comp = calc_av(E,T_comp,p_comp)

#ratio_cur = 10
J_cur = 0
J_err = 100

Jdiff = 1.
Jdiff_min = 0.00001


act_in = np.zeros((ms,ms))
act_out = np.zeros((ms,ms))
act_delta = np.zeros((ms,ms))
for i in range(ms):
	for j in range(ms):	
		if (T[i,j] > 0.00001):
			act_in[i,j]  =  T[i,j] * T[j,i]  
			act_out[i,j] =  T_comp[i,j] * T_comp[j,i] 

outlist = []
ms_red = ms
for i in range(ms):
	if (np.isclose(p_in[i],0.)):
		outlist.append(i)
		ms_red -=1

print("outlist",outlist)

p_min= np.zeros(ms)
T_min = np.zeros((ms,ms))

T_red = np.delete(T,outlist,axis =0)
T_red = np.delete(T_red,outlist,axis =1)

T_min = np.delete(T,outlist,axis =0)
T_min = np.delete(T_min,outlist,axis =1)
# information of T is copied to T_min for later use here. 

T_red, p_red, c_cur = reweight(T_red,ms_red, act_out)


print(outlist)
for i in range(len(outlist)):
	T_red = np.insert(T_red,(outlist[i]),0.,axis=0)
	T_min = np.insert(T_min,(outlist[i]),0.,axis=0)
	T_red = np.insert(T_red,(outlist[i]),0.,axis=1)
	T_min = np.insert(T_min,(outlist[i]),0.,axis=1)


for i in range(len(outlist)):
	p_red = np.insert(p_red,outlist[i],0.)
	p_min = np.insert(p_min,outlist[i],0.)
	

p_cur = p_red

T_cur= T_red


J_cur = calc_av(r,T_cur,p_cur)
J_min = calc_av(r,T_min,p_min)

np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/T/p_"+str(identin)+"_"+str(ident)+".dat",p_cur)
np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/T/pm_"+str(identin)+"_"+str(ident)+".dat",p_min)

####### only for equ systems
Ev,Evec =  analyse_MSM(T_cur)
Enum =3
Evecr = np.zeros((ms,Enum))

np.savetxt(path0+'single_particle/'+str(int(lag))+'/reweight/EV/Evec_'+str(identin)+'_'+str(ident)+'.dat',np.transpose(Evec).real)

np.savetxt(path0+'single_particle/'+str(int(lag))+'/reweight/EV/EV_'+str(identin)+'_'+str(ident)+'.dat',Ev.real)
#############



Sprod_rew = calc_Sprod(T_cur)
Sprod_min = calc_Sprod(T_min)
MFPT_rew = mfpt_markov(T_cur , minpos, rancut, ms,lvl).reshape(lvl*lvl)
#ratio_rew =  np.sum([MFPT_rew[1]/MFPT_rew[3] , MFPT_rew[6]/MFPT_rew[2] , MFPT_rew[5]/MFPT_rew[7]]) / 3.
Eav_rew = calc_av(E, T_cur, p_cur)


#ratio_comp = sum(np.abs([MFPT_comp[1]/MFPT_comp[3] , MFPT_comp[6]/MFPT_comp[2] , MFPT_comp[5]/MFPT_comp[7] ] ))
#print(calc_Caliber(T_cur,p_cur,T,J_comp,r,ms,gamma_cur))


length = 10000

MFPT_hist = np.zeros((lvl*lvl*4,length))
mom_cur = np.zeros((lvl*lvl,3))
mom_curmin = np.zeros((lvl*lvl,3))
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
			MFPT_hist[i*lvl+j+ 3*lvl*lvl] = first_passage(start,end, T_min, length)


			for k in range(1,4):
				mom_cur[i*lvl+j,k-1] = moment_distribution(k,MFPT_hist[i*lvl+j])
				mom_curmin[i*lvl+j,k-1] = moment_distribution(k,MFPT_hist[i*lvl+j+3*lvl*lvl])
				mom_comp[i*lvl+j,k-1] = moment_distribution(k,MFPT_hist[i*lvl+j+lvl*lvl])
				mom_start[i*lvl+j,k-1] = moment_distribution(k,MFPT_hist[i*lvl+j+2*lvl*lvl])

			mom_cur[i*lvl+j,1] = np.sqrt( mom_cur[i*lvl+j,1] - mom_cur[i*lvl+j,0]* mom_cur[i*lvl+j,0] ) #standard deviation
			mom_cur[i*lvl+j,2] = (mom_cur[i*lvl+j,2] - 3.*mom_cur[i*lvl+j,0] * mom_cur[i*lvl+j,1] *mom_cur[i*lvl+j,1]  - mom_cur[i*lvl+j,0]* mom_cur[i*lvl+j,0]* mom_cur[i*lvl+j,0])/ (mom_cur[i*lvl+j,1] * mom_cur[i*lvl+j,1] * mom_cur[i*lvl+j,1] ) #skewness

			mom_comp[i*lvl+j,1] = np.sqrt( mom_comp[i*lvl+j,1] - mom_comp[i*lvl+j,0]* mom_comp[i*lvl+j,0] ) #standard deviation
			mom_comp[i*lvl+j,2] = (mom_comp[i*lvl+j,2] - 3.*mom_comp[i*lvl+j,0] * mom_comp[i*lvl+j,1] *mom_comp[i*lvl+j,1]  - mom_comp[i*lvl+j,0]* mom_comp[i*lvl+j,0]* mom_comp[i*lvl+j,0])/ (mom_comp[i*lvl+j,1] * mom_comp[i*lvl+j,1] * mom_comp[i*lvl+j,1] ) #skewness




for k in range(1,4):
	np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/MFPT/moment"+str(int(k))+"_"+str(identin)+"_"+str(ident)+".dat",mom_cur[:,k-1].reshape((lvl,lvl)))
	np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/MFPT/moment"+str(int(k))+"_"+str(ident)+".dat",mom_comp[:,k-1].reshape((lvl,lvl)))


act_m = np.zeros((ms,ms))
ratio_m = np.zeros((ms,ms))
ratio = np.zeros((ms,ms))
act = np.zeros((ms,ms))
cutoff = 0.0000001
for k in range(ms):
	for l in range(ms):
		if (T_cur[k,l] > cutoff and T_cur[l,k] > cutoff):
			ratio_m[k,l] = ( T_min[k,l] / T_min[l,k] )
			ratio[k,l] = ( T_cur[k,l] / T_cur[l,k] )
			act_m[k,l] = ( T_min[k,l] * T_min[l,k] )
			act[k,l] = ( T_cur[k,l] * T_cur[l,k] )
	
np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/T/ratio_m_"+str(identin)+"_"+str(ident)+".dat",ratio_m)
np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/T/act_"+str(identin)+"_"+str(ident)+".dat",act)
np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/T/ratio_"+str(identin)+"_"+str(ident)+".dat",ratio)
np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/T/act_m_"+str(identin)+"_"+str(ident)+".dat",act_m)
#np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/T/ratio_diffS_"+str(identin)+"_"+str(ident)+".dat",np.abs(sp-ratio_t))


#print(identin, ident)
#for i in [1,2,3,4,6,7,8,9,11,12,13,14]:
#	print(i+1,mom_comp[i,0], mom_cur[i,0],mom_curmin[i,0] )
for i in [1,2,3,5,6,7]:
	print(i+1,mom_comp[i,0], mom_cur[i,0],mom_curmin[i,0] )



print(identin, ident,J_comp,J_cur,J_min, Sprod_comp, Sprod_rew,Sprod_min)


MFPT_cumhist = hist_to_cumhist(MFPT_hist)
#print(identin, ident,gamma_cur, tot_error)
np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/MFPT/hist_"+str(identin)+"_"+str(ident)+".dat", np.transpose(MFPT_hist))
np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/T/T_"+str(identin)+"_"+str(ident)+".dat", T_cur)
np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/T/Tm_"+str(identin)+"_"+str(ident)+".dat", T_min)


	





