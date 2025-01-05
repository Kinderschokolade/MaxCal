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

def reweight(T,r, gamma, ms):
	#W = [ [ np.power( T[i,j] , 2./3.) * np.power(T[j,i], 1./3.) * math.exp(2./3.*gamma*r[i,j] + 1./3. *gamma*r[j,i]) for j in range(ms)] for i in range(ms)]
	W = [ [  T[i,j] * math.exp(gamma*r[i,j]) for j in range(ms)] for i in range(ms)]
	W = np.asarray(W)

	Ev,Evecl = np.linalg.eig(np.transpose(W))
	idx = Ev.argsort()[::-1] 
	Ev = Ev[idx]  #order eigenvectors by size
	Evecl = np.transpose(Evecl)
	Evecl = Evecl[idx,:] 
	Ev,Evecr = np.linalg.eig(W) # algorithm is called double, do this faster!
	idx = Ev.argsort()[::-1] 
	Ev = Ev[idx]  #order eigenvectors by size
	Evecr = np.transpose(Evecr)
	Evecr = Evecr[idx,:] 
	p  = Evecl[0,:] * Evecr[0,:]
	LMm = np.real(-np.log(Evecr[0,:]))
	LMn = 1. - LMm
	LMz = np.real(-np.log(Ev[0]))
	p = np.absolute(np.asarray(p)/sum(p))
	k = [[ Evecr[0,j] / (Evecr[0,i] * Ev[0]) *W[i,j] for j in range(ms)] for i in range(ms)]
	k= np.asarray(k)

	for i in range(ms):
		k[i] = k[i] / sum(k[i,:])

	return np.real(k),np.real(p), LMn, LMm, LMz


def pij(T,q,r,gamma,ms,Sp,LMn,LMm,LMz,w): # with pij to weight S.
	Tnew = np.zeros((ms,ms))
	cutoff = 0.0001
	#print(LMm)
	for i in range(ms):
		for j in range(ms):
			if (q[i,j] > cutoff and q[j,i] > cutoff):
				if (i !=j):
					#print(i,j,w[i,j],Sp[i,j], q[i,j]**(1.-w[i,j]) , q[j,i]**w[i,j] ,-1.+LMm[i]*(-1.+w[i,j]) -LMm[j] *w[i,j] +LMn[j]*(-1.+w[i,j]) - LMn[j] *w[i,j] - LMz +gamma*(1.-w[i,j])*r[i,j]+ gamma *w[i,j]*r[j,i] + w[i,j] * Sp[i,j]  )
					Tnew[i,j] = q[i,j]**(1.-w[i,j]) * q[j,i]**w[i,j] * np.exp(-1.+LMm[i]*(-1.+w[i,j]) -LMm[j] *w[i,j] +LMn[j]*(-1.+w[i,j]) - LMn[j] *w[i,j] - LMz +gamma*(1.-w[i,j])*r[i,j]+ gamma *w[i,j]*r[j,i] + w[i,j] * Sp[i,j]  )
				else:
					Tnew[i,j] = q[i,j] * np.exp(-1.-LMm[i]-LMn[j]- LMz +gamma*r[i,j] )
	
	for k in range(ms):
		Tnew[k,:] = Tnew[k,:] / sum(Tnew[k,:])
	
				
	return Tnew


def pij2(T,q,r,gamma,ms,Sp,LMn,LMm,LMz,a): # with pii pij to weight Sdot
	Tnew = np.zeros((ms,ms))
	cutoff = 0.00001
	#print(LMm)
	w = -a/(1-a)
	for i in range(ms):
		for j in range(ms):
			if (q[i,j] > cutoff and q[j,i] > cutoff):
				if (i > j):
					Tnew[i,j] = q[i,j]**(1.-w[i,j]) * q[j,i]**w[i,j] * np.exp(-1.+LMz+(LMm[i]+LMn[j]-gamma*r[i,j])*(1.-w[i,j]) + (LMm[j]+LMn[i]-gamma*r[j,i])*w[i,j] - w[i,j] * (1-np.exp(Sp[i,j] -Sp[i,j]) - a[i,j] *w[i,j] * np.exp(Sp[i,j])  ))
					print(i,j,Tnew[i,j])
	#				print(i,j,a[i,j], w[i,j], -1.+LMz+(LMm[i]+LMn[j]-gamma*r[i,j])*(1.-w[i,j]) + (LMm[j]+LMn[i]-gamma*r[j,i])*w[i,j] - w[i,j] * (1-np.exp(Sp[i,j] -Sp[i,j]) - a[i,j] *w[i,j] * np.exp(Sp[i,j])  ))
				else:
					Tnew[i,j] = q[i,j] * np.exp(-1.+LMm[j]+LMn[i] + LMz -gamma*r[j,i] - a[j,i] *np.exp(Sp[i,j]) )
	#				print(i,j,-1.+LMm[j]+LMn[i] + LMz -gamma*r[j,i] - a[j,i] *np.exp(Sp[i,j]) )
	
	for k in range(ms):
		Tnew[k,:] = Tnew[k,:] / sum(Tnew[k,:])
	
				
	return Tnew


def aij(T,ms,LMm,LMn,gamma,r,Sp):
	a = np.zeros((ms,ms))
	cutoff = 0.0001
	for i in range(ms):
		for j in range(ms):
			if (i !=j and T[i,j] > cutoff and T[j,i] > cutoff   ):
				a[i,j] = Sp[i,j] -np.log(T[i,j] / T[j,i])  -LMm[i] + LMm[j] -LMn[j] + LMn[i] + gamma * (r[j,i] - r[i,j] )
			else:
				a[i,j] = 0.
	return a


def wij(pi,T,ms):
	w = np.zeros((ms,ms))
	cutoff = 0.0001
	for i in range(ms):
		for j in range(ms):
			if (i !=j and T[i,j] > cutoff and T[j,i] > cutoff   ):
				flow = pi[i] * T[i,j] / (pi[j] *T[j,i])
				w[i,j] = 1./(1.- flow)
			else:
				w[i,j] = 0.
	return w

def wi(S):
	w = np.zeros((ms,ms))
#	cutoff = 0.0001
	for i in range(ms):
		for j in range(ms):
			if (i !=j  ):
				w[i,j] = 1./(1.- np.exp(S[i,j]))
				if(w[i,j] > 50.):
					w[i,j] = 50.
				if(w[i,j] < -50.):
					w[i,j] = -50. # I dont think this is a good solution for Sdot \approx 0


			else: ## actually inf, but it won't be used
				w[i,j] = 0.


	return w


def lagrange_m(T,q,r,gamma,ms,Sp,LMn,LMz,w):
	A = np.zeros((ms,ms))
	b = np.zeros(ms)
	cutoff = 0.0001
	for i in range(ms):
		A[i,i] = -1 - T[i,i] * w[i,i]
		for j in range(ms):
			A[i,j] = -T[i,j] * w[i,j]
			A[i,i] += T[i,j] * w[i,j] 		
			if (q[j,i] > cutoff and q[i,j] > cutoff):	
				b[i] += T[i,j] * w[i,j] * (np.log(q[j,i]/q[i,j]) +LMn[i] - LMn[j] +gamma*(-r[i,j] +r[j,i])  )
			if (j > i):
				b[i] += T[i,j] * w[i,j] *Sp[i,j] 

		
	#print(A)
	#print(b)
	LMm = np.linalg.solve(A,b)

	return LMm
	

def improve_rew(T,q,r,gamma,ms,Sp,LMn,LMm,LMz):
	k =50
	EV,Evec = analyse_MSM(T)
	pi = np.real(Evec[0,:])
	Tnew = T[:,:]
	w = wi(Sp)
	#print("w", w)
	a = aij(q,ms,LMm,LMn,gamma,r,Sp)
	#for i in range(k):
		#print("pi",pi)
#		print("LMm",LMm)
#		print("LMn",LMn)
#		print("LMz",LMz)
#		print("Sp",Sp)
		#w = wij(pi,Tnew,ms)
	#	print(LMm[1], LMm[15],Tnew[0,3],Tnew[10,12])
	Tnew = pij2(T,q,r,gamma,ms,Sp,LMn,LMm,LMz,a)
		#LMm = lagrange_m(Tnew,q,r,gamma,ms,Sp,LMn,LMz,w)	
		#print("Tnew",Tnew)
	
	EV,Evec = analyse_MSM(Tnew)
	pi = np.real(Evec[0,:])
	#print("J",calc_av(r,Tnew,pi))

	return Tnew,pi


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

def define_pathE(ms,cut, ident, identin): 
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

for line in open("/data/isilon/bause/single_particle/param_"+str(ident)+".dat","r"):
	cont = line.split()
	if(cont[0] == 'dT'):
		dT = float(cont[1])
	if(cont[0] == 'microstates'):
		ms = int(cont[1])
	if(cont[0] == 'T0'):
		Temperature = int(cont[1])
	if(cont[0] == 'extf'):
		force = int(cont[1])


potential = np.loadtxt("/data/isilon/bause/single_particle/potential_ms_"+str(ident)+".dat")
path0 = "/data/isilon/bause/"
ms,cut = read_cut(path0,ident)
ms, minpos = read_minpos(path0,ident)
r = define_flow_direction(ms,cut)

#T = np.loadtxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/T/Topt_"+str(identin)+".dat")
#T_comp = np.loadtxt("/data/isilon/bause/single_particle/"+str(lag)+"/T/T_"+str(ident)+".dat")


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


Ev,Evec =  analyse_MSM(T_comp)
p_comp = Evec[0] / sum(Evec[0]) 
temp = np.loadtxt("/data/isilon/bause/single_particle/potential_ms_"+str(ident)+".dat")
E = np.zeros((ms,ms))

for i in range(ms):
	for j in range(ms):
		E[i,j] = (temp[i,1] + temp[j,1]) / 2.

Ng = 500
lvl = len(minpos)
MFPT_ar = np.zeros((Ng, lvl *lvl))
ratio = np.zeros((Ng,3))

J_ar = np.zeros(Ng)

MFPT_comp = mfpt_markov(T_comp , minpos, rancut, ms,lvl).reshape(lvl*lvl)
ratio_comp =  sum( [MFPT_comp[1]/MFPT_comp[3] , MFPT_comp[6]/MFPT_comp[2] , MFPT_comp[5]/MFPT_comp[7]] ) / 3.
J_comp =  calc_av(r, T_comp, p_comp)


Sprod_comp = calc_Sprod(T_comp)
Eav_comp = calc_av(E,T_comp,p_comp)

gamma_start = -0.1
gamma_cur = gamma_start
#ratio_cur = 10
J_cur = 0
J_err = 100

Sp = np.zeros((ms,ms))
for k in range(ms):
	for l in range(ms):
		if (np.abs(k-l) <15):
				Sp[k,l] = ( force / ms * (k-l)+ potential[k,1] - potential[l,1] ) / Temperature
		elif (l>k):
				Sp[k,l] = ( force / ms * (ms-l+k)+ potential[k,1] - potential[l,1] ) / Temperature
		else:
				Sp[k,l] = ( -force / ms * (ms-k+l)+ potential[k,1] - potential[l,1] ) / Temperature


gamma = np.zeros(Ng)

Tnew = np.zeros((ms,ms))
for i in range(Ng):
	gamma[i] =   i*2*abs(gamma_start)/Ng + gamma_start  #random choice here
#	gamma[i] =   i*0.0005 + gamma_start  #random choice here

	Trew, prew,LMn,LMm,LMz = reweight(T, r, gamma[i] , ms)

	J_ar[i] = calc_av(r,Trew,prew)
	print(gamma[i],J_ar[i])
	Trew,prew = improve_rew(Trew,T,r,gamma[i],ms,Sp,LMn,LMm,LMz)

#	Tnew = np.copy(Trew)
#	pnew = np.copy(prew)
	for k in range(ms):
		Trew[k,:] = Trew[k,:] / sum(Trew[k,:])
		
#	EV,Evec = analyse_MSM(Trew)
#	prew = np.real(Evec[0,:])
	
	#Ev, Evec = analyse_MSM(Trew)
	MFPT_ar[i] = mfpt_markov(Trew, minpos, rancut, ms,lvl).reshape(lvl*lvl)
	ratio[i] = [MFPT_ar[i][1]/MFPT_ar[i][3] , MFPT_ar[i][6]/MFPT_ar[i][2] , MFPT_ar[i][5]/MFPT_ar[i][7] ] 
	J_ar[i] = calc_av(r,Trew,prew)
	Sprod = calc_Sprod(Trew)

		
	print(gamma[i], J_ar[i] , J_comp,Sprod, Sprod_comp)#, sum(sum(Tnew[:,:])), sum(np.matmul(pnew,Tnew)))

	if ( np.abs(J_ar[i] - J_comp) < J_err):
		#ratio_cur = sum(np.abs(ratio[i] - ratio_comp))
		#ratcur = sum(np.abs(ratio[i]))
		J_err = np.abs(J_ar[i] - J_comp)
		gamma_cur = gamma[i]
		T_cur = Trew[:,:]
		p_cur = prew[:]
		J_cur = J_ar[i]

		#print( Trew[0,1], T_cur[0,1])
#print(T_cur)
#print(Trew)

Sprod_rew = calc_Sprod(T_cur)
MFPT_rew = mfpt_markov(T_cur , minpos, rancut, ms,lvl).reshape(lvl*lvl)
ratio_rew =  np.sum([MFPT_rew[1]/MFPT_rew[3] , MFPT_rew[6]/MFPT_rew[2] , MFPT_rew[5]/MFPT_rew[7]]) / 3.
Eav_rew = calc_av(E, T_cur, p_cur)

ratio_t = np.zeros((ms,ms))
cutoff = 0.0001
for k in range(ms):
	for l in range(ms):
		if (t_cur[k,l] > cutoff and t_cur[l,k] > cutoff):
			ratio_t[k,l] = np.log( t_cur[k,l] / t_cur[l,k] )
			if (np.abs(k-l) <15):
					sp[k,l] = ( force / ms * (k-l)+ potential[k,1] - potential[l,1] ) / temperature
			elif (l>k):
					sp[k,l] = ( force / ms * (ms-l+k)+ potential[k,1] - potential[l,1] ) / temperature
			else:
					sp[k,l] = ( -force / ms * (ms-k+l)+ potential[k,1] - potential[l,1] ) / temperature

np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/T/ratio_im_"+str(identin)+"_"+str(ident)+".dat",ratio_t)
np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/T/ratio_th_"+str(identin)+"_"+str(ident)+".dat",sp)


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





print(identin, ident,gamma_cur, Sprod_comp, Sprod_rew,  J_comp,J_cur,  mom_comp[1,0],mom_cur[1,0], mom_comp[1,1],mom_cur[1,1], mom_comp[1,2],mom_cur[1,2])

#tot_error= sum(sum(abs(MFPT_hist[i,j] - MFPT_hist[i+lvl*lvl,j]) for i in range(lvl*lvl)) for j in range(length))

MFPT_cumhist = hist_to_cumhist(MFPT_hist)
#print(identin, ident,gamma_cur, tot_error)
np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/MFPT/hist_im_"+str(identin)+"_"+str(ident)+".dat", np.transpose(MFPT_hist))
np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/MFPT/cumhist_im_"+str(identin)+"_"+str(ident)+".dat", np.transpose(MFPT_cumhist))
np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/T/T_im_"+str(identin)+"_"+str(ident)+".dat", T_cur)

	
gamma = np.transpose([gamma])
J_ar = np.transpose([J_ar])

out = np.concatenate((gamma, J_ar),axis = 1)

np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/J/J_im_"+str(identin)+".dat", out)



