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


def equ_alt(c,ms,W,w):
	tup = [ (sum([W[j,k] *c[k]**w[j,k] * c[j]**w[k,j]  for k in range(ms)] )-1.0) for j in range(ms)]
	tup = tuple(tup)
	return tup

def reweight_alt(T, ms, Sprod):  # reweight without detailed balance, exact soution ; wrong!
	w = [ [  1./(1.+ np.exp(Sprod[i,j])) for j in range(ms)] for i in range(ms)]
	w = np.asarray(np.real(w))
	W = [ [  T[i,j]**w[j,i]  *T[j,i]**w[i,j] * math.exp( Sprod[i,j]/2. ) for j in range(ms)] for i in range(ms)]
	W = np.asarray(np.real(W))	

	Omega = np.zeros((ms,ms)) 
	c = np.ones(ms);
	copt = optimize.fsolve(equ_alt,c,args=(ms,W,w,),xtol=1e-12) # should converge for any starting point
	k = [[copt[j]**w[i,j] *copt[i]**w[j,i] *W[i,j] for j in range(ms)] for i in range(ms)]
	k= np.asarray(k)

	for i in range(ms): # minimal correction from numerics, eta = 1.
		norm = sum(k[i,:])
		k[i,:] = k[i,:] / norm

	Ev,Evec =  analyse_MSM(k)
	p = np.real(Evec[0] / sum(Evec[0]) )
	J = calc_av(r,k,p)
	J_err =  np.abs(J - J_comp)

	for i in range(ms):
		for j in range(ms):
			Omega[i,j] = T[i,j] * T[j,i] * np.sqrt(c[i] + c[j]) 

	norm = sum(Omega[:,:])
	Omega = Omega / norm

	return Omega, np.real(k),np.real(p),copt


def reweight(T, ms, Sprod):
	W = [ [  np.sqrt(T[i,j] *T[j,i]) * math.exp( Sprod[i,j]/2. ) for j in range(ms)] for i in range(ms)]
#	W = [ [  T[i,j] * math.exp( (Sprod[i,j]-Sprod_in[i,j])/2. ) for j in range(ms)] for i in range(ms)]
	W = np.asarray(np.real(W))
	Omega = np.zeros((ms,ms))

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
	c = np.abs(np.real(Evecr[0,:]))
	#c = np.zeros(ms);
	#copt = optimize.fsolve(equ,c,args=(ms,W,),xtol=1e-8) # should converge for any starting point
	copt = optimize.least_squares(equ, c,args=(ms,W,), bounds=([0.]*ms, 100.))
	copt= copt.x
#	copt= copt.x
	k = [[copt[j] *copt[i] *W[i,j] for j in range(ms)] for i in range(ms)]
	k= np.asarray(k)
	for i in range(ms): # minimal correction from numerics, eta = 1.
		norm = sum(k[i,:])
		k[i,:] = k[i,:] / norm

	Ev,Evec =  analyse_MSM(k)
	p = np.real(Evec[0] / sum(Evec[0]) )
	J = calc_av(r,k,p)
	J_err =  np.abs(J - J_comp)

	#c = np.ones(ms)
	#W = [ [  np.sqrt(T[i,j] *T[j,i])  for j in range(ms)] for i in range(ms)]
	#W = np.asarray(W)
	#c = optimize.fsolve(equ,c,args=(ms,W,),xtol=1e-12) # should converge for any starting point
	for i in range(ms):
		for j in range(ms):
			Omega[i,j] = T[i,j] * T[j,i] * np.sqrt(c[i] + c[j]) 

#	copt = optimize.fsolve(equ,c,args=(ms,Omega,),xtol=1e-12) # should converge for any starting point
#	kn = [[np.exp(1./2.*(copt[j] +copt[i] ) ) * Omega[i,j] for j in range(ms)] for i in range(ms)]
#	kn = np.asarray(kn)
#	print(k[14,:])
#	print(kn[14,:])
#
##	np.savetxt("T.dat",k)
#	np.savetxt("Tn.dat",kn)

	norm = sum(Omega[:,:])
	Omega = Omega / norm

	
	return Omega, np.real(k),np.real(p),copt

def equ_pij(c,ms,W,w):
	tup = [ (sum([W[i,j] * np.exp((-2.*w[i,j]+1.) * c[i] + (-1.+2.*w[i,j])* c[j]  ) for j in range(ms)] )-1.0) for i in range(ms)]
	tup = tuple(tup)
	return tup

def minimiser_pij(x,ms,W,w,Sprod,Sprod_in):
	global p_glob
	for i in range(ms):
		for j in range(ms):
			p_glob[i,j] = np.exp( x[ms]+ (2.*w[i,j]-1.) * x[j] + x[i] * (-2.*w[i,j]+1.) )  *W[i,j]

	func = sum(  abs(sum(p_glob[i,:] ) -1. ) for i in range(ms))
	return func



def reweight_pij(T, ms, Sprod,Sprod_in): # solution where detailed balance is constraint in Caliber (see notes). 
	w = [ [  1./(1.+ np.exp(Sprod[i,j])) for j in range(ms)] for i in range(ms)]
	w = np.asarray(np.real(w))
	W = [ [  T[i,j]**w[j,i] *T[j,i]**w[i,j] * math.exp( w[i,j]* ( Sprod[i,j] ) )  for j in range(ms)] for i in range(ms)]
	W = np.asarray(np.real(W))


	c = np.ones(ms+1);
	Omega, k,p, copt = reweight(T,ms,Sprod)
	c[0:ms] = copt[:] 
	#c = np.random.rand(ms+1)*5
	
	#copt = optimize.fsolve(equ_pij,c,args=(ms,W,w),xtol=1e-6) # should converge for any starting point

	#out = optimize.minimize(minimiser_pij,c,args=(ms,W,w,Sprod,Sprod_in), method="L-BFGS-B", options={"gtol":1e-02,"ftol":1e-6})
	#copt = out.x
	#print(out.status,out.message)		

	minimizer_kwargs = {"method": "L-BFGS-B", "args": (ms,W,w,Sprod,Sprod_in,), "options" : {"gtol":1e-02,"ftol":1e-6}}
	ret = optimize.basinhopping(minimiser_pij, c, minimizer_kwargs=minimizer_kwargs,niter=500,disp="True",stepsize=0.2)
	print(ret)
	copt = ret.x
	k = [[np.exp( copt[ms]+ (2.*w[i,j]-1.) * copt[j] + copt[i] * (-2.*w[i,j]+1.) )  *W[i,j] for j in range(ms)] for i in range(ms)]
	k= np.asarray(k)

	norm= np.ones(ms+1)
	for i in range(ms): # minimal correction from numerics, eta = 1.
		norm[i] = sum(k[i,:])
		k[i,0:ms] = k[i,:] / norm[i]
	#np.savetxt("c.dat",np.transpose( np.concatenate(([copt],[norm]),axis =0 )))

	Ev,Evec =  analyse_MSM(k)
	p = np.real(Evec[0] / sum(Evec[0]) )

	Omega = np.zeros((ms,ms))
	return Omega, np.real(k),np.real(p),copt




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

def minimiser_1st(x,ms,W,Sprod,Sprod_in):
	global p_glob
	global alpha_glob
	for i in range(ms):
		alpha_glob[i] = 0.
		for j in range(ms):
			p_glob[i,j] =W[i,j]*np.exp( Sprod[i,j]/4. *( x[ms+i] - x[ms+j] )  + 1./2. * (x[i] + x[j]) )
			alpha_glob[i] += p_glob[i,j] * ( 1./2. *(Sprod[i,j] - Sprod_in[i,j] - x[ms+i] + x[ms+j])  + Sprod[i,j]/4.*( x[ms+i] - x[ms+j] - Sprod_in[i,j]) )

	func1 =  sum( abs(sum(p_glob[i,:]) -1. ) for i in range(ms)) #+ 0.0001 *np.linalg.norm(x)
	func2 =  sum( abs(x[i] - x[2*ms] + alpha_glob[i]-1.) for i in range(ms)) #+ 0.0001 *np.linalg.norm(x)
	func = func1+func2
	return func

def minimiser_2nd(x,ms,W,Sprod,Sprod_in):
	global p_glob
	global alpha_glob
	for i in range(ms):
		alpha_glob[i] = 0.
		for j in range(ms):
			p_glob[i,j] =W[i,j]*np.exp( -Sprod[i,j]/4. *( -x[ms+i] + x[ms+j] +Sprod[i,j] - Sprod_in[i,j] ) + 1./2. * (x[i] + x[j]) )
			alpha_glob[i] += p_glob[i,j] * ( 1./2. *(Sprod[i,j] - Sprod_in[i,j] - x[ms+i] + x[ms+j])  - Sprod[i,j]/4.*( -x[ms+i] + x[ms+j] +Sprod[i,j] - Sprod_in[i,j]) )

	func1 =  sum( abs(sum(p_glob[i,:]) -1. ) for i in range(ms)) #+ 0.0001 *np.linalg.norm(x)
	func2 =  sum( abs(x[i] - x[2*ms] + alpha_glob[i]-1.) for i in range(ms)) #+ 0.0001 *np.linalg.norm(x)
	func = func1+func2
	return func


def minimiser_3rd(x,ms,W,Sprod,Sprod_in):
	global p_glob
	global alpha_glob
	for i in range(ms):
		alpha_glob[i] = 0.
		for j in range(ms):
			p_glob[i,j] =W[i,j]*np.exp( -Sprod[i,j]/4. *( -x[ms+i] + x[ms+j] ) + Sprod[i,j]**(3) /48. *(x[ms+i] -x[ms+j])   + 1./2. * (x[i] + x[j]) )
			alpha_glob[i] += p_glob[i,j] * ( 1./2. *(Sprod[i,j] - Sprod_in[i,j] - x[ms+i] + x[ms+j])  - Sprod[i,j]/4.*( -x[ms+i] + x[ms+j]) +Sprod[i,j]**(3)/48. *(x[ms+j] - x[ms+i]) )

	func1 =  sum( abs(sum(p_glob[i,:]) -1. ) for i in range(ms)) #+ 0.0001 *np.linalg.norm(x)
	func2 =  sum( abs(x[i] - x[2*ms] + alpha_glob[i]-1.) for i in range(ms)) #+ 0.0001 *np.linalg.norm(x)
	func = func1+func2
	return func

def minimiser_3rd_full(x,ms,W,Sprod,Sprod_in):
	global p_glob
	global alpha_glob
	for i in range(ms):
		alpha_glob[i] = 0.
		for j in range(ms):
			p_glob[i,j] =W[i,j]*np.exp( -Sprod[i,j]/4. *( -x[ms+i] + x[ms+j] +Sprod[i,j] - Sprod_in[i,j] ) + Sprod[i,j]**(3) /48. *(x[ms+i] -x[ms+j] -Sprod_in[i,j])   + 1./2. * (x[i] + x[j]) )
			alpha_glob[i] += p_glob[i,j] * ( 1./2. *(Sprod[i,j] - Sprod_in[i,j] - x[ms+i] + x[ms+j])  - Sprod[i,j]/4.*( -x[ms+i] + x[ms+j] +Sprod[i,j] - Sprod_in[i,j]) +Sprod[i,j]**(3)/48. *(x[ms+j] - x[ms+i] -Sprod_in[i,j]) )

	func1 =  sum( abs(sum(p_glob[i,:]) -1. ) for i in range(ms)) #+ 0.0001 *np.linalg.norm(x)
	func2 =  sum( abs(x[i] - x[2*ms] + alpha_glob[i]-1.) for i in range(ms)) #+ 0.0001 *np.linalg.norm(x)
	func = func1+func2
	return func


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

def minimiser_trafo2(x,ms,W,w,Sprod,Sprod_in):
	global p_glob
	global alpha_glob
	for i in range(ms):
		alpha_glob[i] = 0.
		for j in range(ms):
			term[i,j] =  ( Sprod[i,j]-Sprod_in[i,j] + x[ms+j] - x[ms+i])
			p_glob[i,j] =W[i,j]*np.exp( (w[i,j]-w[j,i])/2.*term[i,j]  + 1./2. * (x[i] + x[j]) )
			alpha_glob[i] += p_glob[i,j] *w[i,j] * term[i,j]

	func1 =  sum( np.abs(sum(p_glob[i,:]) -1. ) for i in range(ms)) #+ 0.0001 *np.linalg.norm(x)
	func2 =  sum( np.abs(x[i] - x[2*ms] + alpha_glob[i]-1.) for i in range(ms)) #+ 0.0001 *np.linalg.norm(x)
	func = 2.*func1+func2
	return func

def reweight_full2(T,ms,Sprod):
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
			W[i,j] = T[i,j] *np.exp( 1./2. *(Sprod[i,j] - Sprod_in[i,j] ) - 1.)	

	x = np.zeros(2*ms+1)
	x[:ms] = copt[:]

	np.savetxt("c.dat", copt)
	for ii in range(1):
		minimizer_kwargs = { "args": (ms,W,w,Sprod,Sprod_in),"method":"L-BFGS-B", "options":{"gtol": 1e-2,"ftol":1e-6}, }
		#out = optimize.basinhopping(minimiser_trafo2, x,minimizer_kwargs=minimizer_kwargs, callback = callback, disp = 'True', niter = 20, interval = 100 , stepsize=1.0 )
		out = optimize.minimize(minimiser_trafo2,x,args=(ms,W,w,Sprod,Sprod_in), method="L-BFGS-B", options={"gtol":1e-02,"ftol":1e-6})
		x = out.x
		v = x[:ms]
		u = x[ms:2*ms]
		zeta =x[2*ms]
		k = [ [ W[i,j]*np.exp( (w[i,j]-w[j,i] )/2. * (Sprod[i,j] - Sprod_in[i,j] +u[j] - u[i])+ 1./2.* (v[i] + v[j] )  ) for j in range(ms) ] for i in range(ms)]	
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
				W[i,j] = T[i,j] *np.exp( 1./2. *(Sprod[i,j] - Sprod_in[i,j] ) - 1.)	

	np.savetxt("u.dat",u)
	np.savetxt("v.dat",v)

	return Invariant,np.real(k),np.real(p), x


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




def reweight_minimise(T,ms,Sprod):
	Omega,k,p,copt = reweight(T, ms, Sprod)
	w = np.zeros((ms,ms))
	cutoff = 0.0001
	for i in range(ms):
		for j in range(ms):
			if (k[j,i] > cutoff):
				w[i,j] = 1./ (1+ p[i] / p[j] * np.exp(Sprod[i,j]) )
			else:
				w[i,j] = 1./2.	

	W = np.zeros((ms,ms))
	W_approx = np.zeros((ms,ms))
	Sprod_in = np.zeros((ms,ms))
	for i in range(ms):
		for j in range(ms):
#			W[i,j] = T[i,j]**w[j,i] * T[j,i]**w[i,j] *np.exp(w[i,j]* Sprod[i,j] -1.)	
			if(T[i,j] > cutoff and T[j,i] > cutoff):
				Sprod_in[i,j] = np.log( T[i,j] / T[j,i] )
			W[i,j] = T[i,j] *np.exp( w[i,j] *(Sprod[i,j] - Sprod_in[i,j] ) - 1.)	
			W_approx[i,j] = T[i,j] *np.exp( 1./2. *(Sprod[i,j] - Sprod_in[i,j]) - 1.)	

	xhist = np.zeros((2*ms+1,6))
	Omega,k,p,copt = reweight(T, ms, Sprod)
	x = np.zeros(2*ms+1)
	x[:ms] = copt[:]
	xhist[:ms,0] = copt[:]
	minimizer_kwargs = { "args": (ms,W,w,Sprod,Sprod_in),"method":"L-BFGS-B",}# "options":{"gtol": 1e-4,"ftol":1e-6}, }
#	out = optimize.basinhopping(minimiser_1st, x,minimizer_kwargs=minimizer_kwargs, callback = callback, disp = 'True', niter = 10, interval = 100 , stepsize=0.5 )

	act = np.zeros((ms,ms))

	out = optimize.minimize(minimiser_1st,x,args=(ms,W_approx,Sprod,Sprod_in), method="L-BFGS-B")
	print("optimizer 1st:",out.fun)
	x = out.x
	v = x[:ms]
	u = x[ms:2*ms]
	zeta =x[2*ms]
#	xhist[:,1] = out.x
	k = [ [ W_approx[i,j]*np.exp(Sprod[i,j]/4. *(-u[i]+u[j] )  +1./2.* (v[i] + v[j] )  ) for j in range(ms) ] for i in range(ms)]	
	act = np.asarray([[ np.log( k[i][j] * k[j][i]) for j in range(ms)] for i in range(ms)])
	np.savetxt("act1.dat",act)

	k = np.asarray(k)
	Ev,Evec =  analyse_MSM(k)
	p = np.real(Evec[0] / sum(Evec[0]) )	
	for i in range(ms):
		for j in range(ms):
			if (k[j,i] > cutoff):
				w[i,j] = 1./ (1+ p[i] * k[i,j] / (p[j]*k[j,i]))
			else:
				w[i,j] = 1./2.	


	out = optimize.minimize(minimiser_2nd,x,args=(ms,W_approx,Sprod,Sprod_in), method="L-BFGS-B")
	print("optimizer 2nd:",out.fun)
	x = out.x
	v = x[:ms]
	u = x[ms:2*ms]
	zeta =x[2*ms]	
	xhist[:,2] = out.x
	k = [ [ W_approx[i,j]*np.exp( - Sprod[i,j]/4. *(-u[i]+u[j] + Sprod[i][j] - Sprod_in[i,j] )  +1./2.* (v[i] + v[j] )  ) for j in range(ms) ] for i in range(ms)]	
	act = np.asarray([[ np.log( k[i][j] * k[j][i] ) for j in range(ms)] for i in range(ms)])
	np.savetxt("act2.dat",act)
	k = np.asarray(k)
	Ev,Evec =  analyse_MSM(k)
	p = np.real(Evec[0] / sum(Evec[0]) )	
	for i in range(ms):
		for j in range(ms):
			if (k[j,i] > cutoff):
				w[i,j] = 1./ (1+ p[i] * k[i,j] / (p[j]*k[j,i]))
			else:
				w[i,j] = 1./2.	




	out = optimize.minimize(minimiser_3rd,x,args=(ms,W_approx,Sprod,Sprod_in), method="L-BFGS-B")
	print("optimizer 3rd:",out.fun)#
	x = out.x
	v = x[:ms]
	u = x[ms:2*ms]
	zeta =x[2*ms]	
	xhist[:,3] = out.x
	k = [ [ W_approx[i,j]*np.exp( - Sprod[i,j]/4. *(-u[i]+u[j]) +Sprod[i,j]**(3)/48. *(-u[i] +u[j])  +1./2.* (v[i] + v[j] )  ) for j in range(ms) ] for i in range(ms)]	
	act = np.asarray([[ np.log( k[i][j] * k[j][i] ) for j in range(ms)] for i in range(ms)])
	np.savetxt("act3.dat",act)
	k = np.asarray(k)
	Ev,Evec =  analyse_MSM(k)
	p = np.real(Evec[0] / sum(Evec[0]) )	
	for i in range(ms):
		for j in range(ms):
			if (k[j,i] > cutoff):
				w[i,j] = 1./ (1+ p[i] * k[i,j] / (p[j]*k[j,i]))
			else:
				w[i,j] = 1./2.	



#	out = optimize.minimize(minimiser_3rd_full,x,args=(ms,W_approx,Sprod,Sprod_in), method="L-BFGS-B")
#	print("optimizer 3rd full:",out.fun)
#	x = out.x
#	xhist[:,4] = out.x

	out = optimize.minimize(minimiser_trafo,x,args=(ms,W,w,Sprod,Sprod_in), method="L-BFGS-B")
	print("optimizer:",out.fun)

	xhist[:,4] = out.x

	x = out.x

	np.savetxt("xhist.dat", xhist)
	#x = np.asarray(out.x)
	#k = [ [ W[i,j]*np.exp(w[j,i] *(m[i]+n[j] )  +w[i,j] * (m[j] + n[i] )  ) for j in range(ms) ] for i in range(ms)]	
	v = x[:ms]
	u = x[ms:2*ms]
	zeta =x[2*ms]
#	print(zeta)
	#full
	k = [ [ W[i,j]*np.exp(w[i,j] *(-u[i]+u[j] )  +1./2.* (v[i] + v[j] - u[j] + u[i] )  ) for j in range(ms) ] for i in range(ms)]	
	#1st
#	k = [ [ W_approx[i,j]*np.exp(Sprod[i,j]/4. *(-u[i]+u[j] )  +1./2.* (v[i] + v[j] )  ) for j in range(ms) ] for i in range(ms)]	
	#2nd
#	k = [ [ W_approx[i,j]*np.exp( - Sprod[i,j]/4. *(-u[i]+u[j] + Sprod[i][j] - Sprod_in[i,j] )  +1./2.* (v[i] + v[j] )  ) for j in range(ms) ] for i in range(ms)]	
	#3rd
	#k = [ [ W_approx[i,j]*np.exp( - Sprod[i,j]/4. *(-u[i]+u[j]) +Sprod[i,j]**(3)/48. *(-u[i] +u[j])  +1./2.* (v[i] + v[j] )  ) for j in range(ms) ] for i in range(ms)]	
	k = np.asarray(k)
	knorm = np.zeros(ms)
	
	print("norm",[sum( k[i,:]) for i in range(ms)]  ) # should be 1.
	print("glob",[ -zeta+ +v[i] +sum( w[i,j] * k[i,j]*(Sprod[i,j] - Sprod_in[i,j]  -u[i] + u[j] )  for j in range(ms)) for i in range(ms) ] ) # should be 1

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
	p = np.real(Evec[0] / sum(Evec[0]) )
	J = calc_av(r,k,p)


	return Invariant,np.real(k),np.real(p),x

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
path = "/data/isilon/bause/single_particle/"+str(lag)+"/EV/ps_"+str(identin) +".dat"
p_in = np.loadtxt(path)
ms = np.shape(T)[0]
minT = np.nanmin(T)
for i in range(ms):
	for j in range(ms):
		if (T[i,j]  > 0. ):
			T[i,j] = T[i,j] 
		else:
			T[i,j] = minT
		#	T[i,j] = 0. # doesnt make a difference! minT might help to find data wher pij = 0


path0 = "/data/isilon/bause/"
ms,cut = read_cut(path0,ident)
ms, minpos = read_minpos(path0,ident)

#------------------------
minpos = np.asarray([1,9,16,23])
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

sp = np.zeros((ms,ms))
sp_in = np.zeros((ms,ms))

A = np.zeros((ms,ms))

sp_loc = np.zeros(ms)
#----------------------------
f_loc = np.zeros(ms)
#f_loc[10:16] = +force /ms  # driving from 0.312 to 0.5
f_loc[0:ms] = force/ms

#f_loc[19] = 5.739e-05 
#f_loc[20] = 0.0041148 
#f_loc[21] = 0.0870384 
#f_loc[22] = 0.54316 
#f_loc[23] = 1.  
#f_loc[24] = 0.54316 
#f_loc[25] = 0.0870384
#f_loc[26] = 0.0041148 
#f_loc[27] = 5.73909e-05 
#f_loc[28] = 2.36152e-07 


#f_loc[1]=  1.90367e-05  
#f_loc[2]=  0.000535806  
#f_loc[3]=  0.00819126  
#f_loc[4]=  0.0680176  
#f_loc[5]=  0.306775  
#f_loc[6]=  0.75153  
#f_loc[7]=  1  
#f_loc[8]=  0.722739 
#f_loc[9]=  0.283721 
#f_loc[10]=  0.0604961 
#f_loc[11]=  0.00700636 
#f_loc[12]=  0.000440742 
#f_loc[13]=  1.50593e-05 
#
#
#f_loc[:] = f_loc[:] * force/ms

#---------------------------

for i in range(ms-1):
	sp_loc[i] = (potential[i,1] - potential[i+1,1] + f_loc[i] ) /temperature

sp_loc[ms-1] = (potential[ms-1,1] - potential[0,1] + f_loc[ms-1] ) /temperature

for k in range(ms):
	for l in range(ms):
		if (np.abs(k-l) < (ms/2)):
			#sp[k,l] = ( -force / ms * (k-l) + potential[k,1] - potential[l,1] ) / temperature
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
			sp[k,l] = sum(sp_loc[0:l]) +sum(sp_loc[k:ms])
			#sp[k,l] = ( force / ms * (ms-k+l)*mul+ potential[k,1] - potential[l,1] ) / temperature
			#sp_in[k,l] = (force_in / ms * (ms-k+l)*mul+ potential_in[k,1] - potential_in[l,1] ) / temperature_in


#gamma_cur = fmin(reweight_min,x0=0.0,ftol = Jdiff_min , disp=False, args=(T ,r,E , ms , J_comp, sp, sp_in))
#gamma_cur = gamma_cur[0]


#for i in range(3):
#	gamma_cur = i * -0.01
#	J_err,T_cur,p_cur = reweight(gamma_cur,T ,r , ms , J_comp, sp)
#	print(gamma_cur,J_err)

gamma_cur =0

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
sp_red = np.delete(sp,outlist,axis =0)
T_red = np.delete(T_red,outlist,axis =1)
sp_red = np.delete(sp_red,outlist,axis =1)

T_min = np.delete(T,outlist,axis =0)
sp_in_red = np.delete(sp_in,outlist,axis =0)
T_min = np.delete(T_min,outlist,axis =1)
sp_in_red = np.delete(sp_in_red,outlist,axis =1)
# information of T is copied to T_min for later use here. 

#Omega,T_cur, p_cur = reweight_selfit(T,ms,sp,sp_in,delta_min) # does not converge

#Omega,T_min, p_min,c_cur = reweight_minimise(T_min,ms_red,sp_red)
#Omega,T_min, p_min,c_cur = reweight_full2(T_min,ms_red,sp_red)

Omega,T_red, p_red,c_cur = reweight(T_red,ms_red,sp_red)

#Omega,T_min, p_min,c_cur = reweight_A(T_min,ms_red,A)

#Omega,T_min, p_min,c_cur = reweight_pij(T_min,ms_red,sp_red,sp_in_red) # still not as good
#for i in range(1,10):
gamma = 0.

#Omega,T_min, p_min = reweight_minimise(T,ms,sp,gamma)  


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

np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/T/Invariant_"+str(identin)+"_"+str(ident)+".dat",Omega)
np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/T/Invariantn_"+str(identin)+"_"+str(ident)+".dat",E/sum(E))


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
#	np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/MFPT/moment"+str(int(k))+"_"+str(ident)+".dat",mom_comp[:,k-1].reshape((lvl,lvl)))

#
#with open("/data/isilon/bause/single_particle/"+str(lag)+"/Sprod_compo_"+str(identin)+".dat", 'ab') as f:
##with open("/data/isilon/bause/single_particle/"+str(lag)+"/Sprod_comp.dat", 'ab') as f:
#	out = np.zeros(4)
#	out[0] = identin
#	out[1] = ident
#	out[2] = Sprod_comp
#	out[3] = Sprod_rew
#	np.savetxt(f,[out])
#
#with open("/data/isilon/bause/single_particle/"+str(lag)+"/moment1_compo_"+str(identin)+".dat", 'ab') as f:
##with open("/data/isilon/bause/single_particle/"+str(lag)+"/moment1_comp.dat", 'ab') as f:
#	out = np.zeros(2+lvl*lvl*2)
#	out[0] = identin
#	out[1] = ident
#	out[2:lvl*lvl+2] = mom_comp[:,0]
#	out[lvl*lvl+2:] = mom_cur[:,0]
#	np.savetxt(f,[out])
#
#
#with open("/data/isilon/bause/single_particle/"+str(lag)+"/moment2_compo_"+str(identin)+".dat", 'ab') as f:
##with open("/data/isilon/bause/single_particle/"+str(lag)+"/moment2_comp.dat", 'ab') as f:
#	out = np.zeros(2+lvl*lvl*2)
#	out[0] = identin
#	out[1] = ident
#	out[2:lvl*lvl+2] = mom_comp[:,1]
#	out[lvl*lvl+2:] = mom_cur[:,1]
#	np.savetxt(f,[out])
#
#with open("/data/isilon/bause/single_particle/"+str(lag)+"/moment3_compo_"+str(identin)+".dat", 'ab') as f:
##with open("/data/isilon/bause/single_particle/"+str(lag)+"/moment3_comp.dat", 'ab') as f:
#	out = np.zeros(2+lvl*lvl*2)
#	out[0] = identin
#	out[1] = ident
#	out[2:lvl*lvl+2] = mom_comp[:,2]
#	out[lvl*lvl+2:] = mom_cur[:,2]
#	np.savetxt(f,[out])
#

#f = open("/data/isilon/bause/single_particle/"+str(lag)+"/MFPT1_comp.dat", 'a')
#np.savetxt(f,np.array([int(identin),int(ident),mom_comp[1],momrew[1]]))

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


	





