import numpy as np
import math
import random
from scipy.optimize import fsolve
from fct_markov import analyse_MSM

def reweight_J(T,r, gamma, ms):
	#W = [ [ np.power( T[i,j] , 2./3.) * np.power(T[j,i], 1./3.) * math.exp(2./3.*gamma*r[i,j] + 1./3. *gamma*r[j,i]) for j in range(ms)] for i in range(ms)]
	W = [ [  T[i,j] * math.exp(gamma*r[i,j]) for j in range(ms)] for i in range(ms)]
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


def equ(c,ms,W):
	tup = [ (sum([W[j,k] *c[k] * c[j] for k in range(ms)] )-1.0) for j in range(ms)]
	tup = tuple(tup)
	return tup

def reweight(T, ms, Sprod):
	W = [ [  np.sqrt(T[i,j] *T[j,i]) * math.exp( Sprod[i,j]/2. ) for j in range(ms)] for i in range(ms)]
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
	c = np.ones(ms);
	copt = fsolve(equ,c,args=(ms,W,),xtol=1e-12) # should converge for any starting point
	k = [[copt[j] *copt[i] *W[i,j] for j in range(ms)] for i in range(ms)]
	k= np.asarray(k)

	for i in range(ms): # minimal correction from numerics, eta = 1.
		norm = sum(k[i,:])
		k[i,:] = k[i,:] / norm

	Ev,Evec =  analyse_MSM(k)
	p = np.real(Evec[0] / sum(Evec[0]) )

	for i in range(ms):
		for j in range(ms):
			Omega[i,j] = T[i,j] * T[j,i] * np.sqrt(c[i] + c[j]) 

	norm = sum(Omega[:,:])
	Omega = Omega / norm

	
	return Omega, np.real(k),np.real(p)


ms = 4
q = np.zeros((ms,ms))
E = np.zeros(ms)
A = np.zeros((ms,ms))

E[0] =  0.01; E[1] = 0.2 ; E[2] =0.18 ; E[3] =0.02 #;E[4] =0.19; E[5] = 0.05 # Energy levels
stat_p = np.zeros(ms)
beta =1

A[0][1] = 0.2 ; A[0][2] = 0.0 ; A[0][3] = 0.1 ; #A[0][4] = 0.03 ; A[0][5] = 0.09 ;
#A[0][1] = 0.1 ; A[0][2] = 0.0 ; A[0][3] = 0.1 ; #A[0][4] = 0.03 ; A[0][5] = 0.09 ;
A[1][2] = 0.05 ; A[1][3] = 0.0 ; #A[1][4] = 0. ; A[1][5] = 0.05 ;
#A[1][2] = 0.1 ; A[1][3] = 0.0 ; #A[1][4] = 0. ; A[1][5] = 0.05 ;
A[2][3] = 0.1 ;# A[2][4] = 0.06 ; A[2][5] = 0. ;
#A[2][3] = 0.1 ;# A[2][4] = 0.06 ; A[2][5] = 0. ;
#A[3][4] = 0.08 ; A[3][5] = 0.05 ;
#A[4][5] = 0.13
Sp = np.zeros((ms,ms))
for k in range(ms):
	for l in range(ms):
			if (np.abs(k-l) <( ms /2)):
				Sp[k,l] = ( E[k] - E[l] ) *beta
			elif (l>k):
				Sp[k,l] = ( E[k] - E[l] ) *beta
			else:
				Sp[k,l] = ( E[k] - E[l] ) *beta


for i in range(ms):
	stat_p[i] = np.exp(-beta*E[i])

	for j in range(i):
		A[i,j] = A[j,i] 
		q[j,i] = A[i,j] * np.exp(-beta/2.*(E[i] - E[j]))
		q[i,j] = A[i,j] * np.exp(-beta/2.*(E[j] - E[i]))
stat_p = stat_p / sum(stat_p[:])
print("Boltzmann", stat_p)
for i in range(ms):
	q[i][i] = 1. - sum(q[i][:])

if (np.any(q < 0)):
	print(q)

Ev,Evecl = np.linalg.eig(np.transpose(q))
idx = Ev.argsort()[::-1] 
Ev = Ev[idx]  #order eigenvectors by size
Evecl = np.transpose(Evecl)
Evecl = Evecl[idx,:] 

ps = np.real(Evecl[0,:] / sum(Evecl[0,:]))
#print(stat_p,sum(stat_p)) 


pfull = np.zeros((100,ms))
pfulla = np.zeros((100,ms))
pJall = np.zeros((100,ms))

Sfull = np.zeros((100,ms))
Afull = np.zeros((100,ms))

pfull[0] = ps;
pfulla[0] = ps;
pJall[0] = ps;
palt= np.copy(q)
pJ= np.copy(q)
r = np.zeros((ms,ms))

r[0,1] = -1 ; r[0,2] = -2 ; r[0,3] = 1;
r[1,2] = -1; r[1,3] = -2; r[1,0] = -r[0,1]
r[2,3] = -1; r[2,0] = -r[0,2]; r[2,1] = -r[1,2];
r[3,0] = -1; r[3,1] = +2 ; r[3,2] = 1;


for f in range(1,100):
	force = f *0.5
	for k in range(ms):
		for l in range(ms):
			if (np.abs(k-l) <( ms /2)):
				Sp[k,l] = ( -force / ms * (k-l) + E[k] - E[l] ) *beta
				palt[k,l] = A[k,l] * np.exp(beta/2.*(E[k] - E[l] - force* (k-l) /ms ))
			elif (l>k):
				Sp[k,l] = ( -force / ms * (ms-l+k)+ E[k] - E[l] ) *beta
				palt[k,l] = A[k,l] * np.exp(beta/2.*(E[k] - E[l] - force* (ms-l+k) /ms ))
			else:
				Sp[k,l] = ( +force / ms * (ms-k+l)+ E[k] - E[l] ) *beta
				palt[k,l] = A[k,l] * np.exp(beta/2.*(E[k] - E[l] + force* (ms-k+l) /ms ))

	for k in range(ms):
		palt[k,k] =0.	
		pJ[k,k] =0.	
		palt[k][k] = 1. - sum(palt[k][:])
		pJ[k][k] = 1. - sum(pJ[k][:])


	Ev,Evecl = np.linalg.eig(np.transpose(palt))
	idx = Ev.argsort()[::-1] 
	Ev = Ev[idx]  #order eigenvectors by size
	Evecl = np.transpose(Evecl)
	Evecl = Evecl[idx,:] 

	psalt = np.real(Evecl[0,:] / sum(Evecl[0,:]))
#
	Inv,p,ps = reweight(q,ms,Sp)
	gamma = f*0.05
	pJ,psJ = reweight_J(q,r,gamma,ms)

	Sfull[f][0] = Sp[0,1]
	Afull[f][0] = p[0,1] * p[1,0]
	Sfull[f][1] = Sp[1,2]
	Afull[f][1] = p[1,2] * p[2,1]
	Sfull[f][2] = Sp[2,3]
	Afull[f][2] = p[2,3] * p[3,2]
	Sfull[f][3] = Sp[3,0]
	Afull[f][3] = p[3,0] * p[0,3]


	pJall[f] = psJ
	pfull[f] = ps
	#print(palt[0,0])
	if (np.any(palt < 0)):
		#print(f,palt)
		pass
	else:
		pfulla[f] = psalt

#	Teff =  (E[2] - E[5] ) / (np.log(ps[5] / ps[2])) 

acc = np.zeros((ms,ms))

for i in range(ms):
	for j in range(ms):
		A[i,j] = p[i,j] * p[j,i]
		acc[i,j] = np.sqrt( A[i,j] * np.exp(Sp[i,j]) )

#	acc[i,:] /= sum(acc[i,:])
	
# MC simulation at highest force
psim = np.zeros(ms)
Nsim = 10
pos =0;
hist = np.zeros(ms)
T =  np.zeros((ms,ms))

for i in range(Nsim):
	x = random.random()
	pos_new = random.randint(0,3)
	if (x < acc[pos,pos_new]):#A[pos,pos_new]* Sp[pos,pos_new]):
		pos = pos_new
		hist[pos]+=1	
	else:
		hist[pos] +=1

print("MC", hist/sum(hist))
print("rew",ps)	


np.savetxt("pfull.dat", pfull)
np.savetxt("Afull.dat", Afull)
np.savetxt("Sfull.dat", Sfull)
np.savetxt("palt.dat", pfulla)
np.savetxt("pJ.dat", pJall)
 









