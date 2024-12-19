import numpy as np
import argparse
import getpass
import os.path
from define_flow import read_minpos
from define_flow import read_cut
from fpt_ana import first_passage
from msmtools.analysis import mfpt
from define_flow import define_flow_direction
import math

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



def analyse_MSM(k):	
	Ev,Evec = np.linalg.eig(np.transpose(k))
	idx = Ev.argsort()[::-1] 
	Ev = Ev[idx]  #order eigenvectors by size
	Evec = np.transpose(Evec)
	Evec = Evec[idx,:] 
	p = Evec[0] / sum(Evec[0])
	return Ev,p,Evec



def calc_av(F, T, p):
	ms = np.shape(T)[0]
	return np.real( sum( sum( p[i] * T[i,j] *F[i,j] for j in range(ms)) for i in range(ms)  ))


def equ(c,ms,W):
	tup = [ (sum([W[j,k] * np.exp(1./2.*(c[k] + c[j])) for k in range(ms)] )-1.0) for j in range(ms)]
	tup = tuple(tup)
	return tup

def reweight(T, ms, Sprod):
#	W = [ [  T[i,j] * math.exp(gamma/2.*(r[i,j]+r[j,i]) + (Sprod[i,j]-Sprod_in[i,j])/2. ) for j in range(ms)] for i in range(ms)]
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

	copt = optimize.fsolve(equ,c,args=(ms,W,),xtol=1e-12) # should converge for any starting point
#	copt = optimize.least_squares(equ, c,args=(ms,W,), bounds=([0.]*ms, 1.))
#	copt= copt.x
	k = [[np.exp(1./2.*(copt[j] +copt[i] ) ) *W[i,j] for j in range(ms)] for i in range(ms)]
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
			Omega[i,j] = np.sqrt( T[i,j] * T[j,i] ) * np.exp((copt[i] + copt[j])/4.) # or something like that
	
	#print("glob",[ -copt[i] - sum( 1./2. * k[i,j]*(Sprod[i,j] - Sprod_in[i,j] -u[i] + u[j])  for j in range(ms)) for i in range(ms) ] )
	
	norm = sum(Omega[:,:])
	Omega = Omega / norm

	return Omega, np.real(k),np.real(p),copt





parser = argparse.ArgumentParser()
parser.add_argument('-i', help='ident start')
parser.add_argument('-cstart', help='start to count x (kkxk)')
parser.add_argument('-cstop', help='stop to count x (kkxk)')
parser.add_argument('-bstart', help='stop to count x (kxkk)')
parser.add_argument('-bstop', help='stop to count x (kxkk)')
#parser.add_argument('-d', help='number of sets, count up x (kkkx)')
parser.add_argument('-l', help='lagtime chosen')
parser.add_argument('-o', help='reweight all data here')
parser.add_argument('-rc', help='rancut')
args = parser.parse_args();
ident = int(args.i)

cmin = int(args.cstart)
cmax = int(args.cstop)
count_full = cmax- cmin

bmin = int(args.bstart)
bmax = int(args.bstop)
bcount = bmax - bmin

l = int(args.l)
lstr = str(args.l)
aim = str(args.o)
aimc = aim
aimc= aim[0:2]+"0"+aim[3]
rancut = int(args.rc)


user = getpass.getuser()

#if (user=="mbause"):
path0 = "/data/isilon/bause/"
path1 = "/data/pckr194/bause/"


ms, minpos = read_minpos(path0,ident)
ms, cut = read_cut(path0,ident)
lvl = len(minpos)

F = open(path0+"single_particle/count_"+ str(ident) +".dat","r") 
line = F.readline()
lagtime = line.split()
lagtime.pop(0)
lagtime  = [int(i) for i in lagtime]

cols = len(lagtime)
for i in range(cols):
	if (lagtime[i] == l):
		col = i

potential = np.loadtxt("/data/isilon/bause/single_particle/potential_ms_"+str(ident)+".dat")

if (os.path.isfile(path0+"single_particle/param_"+str(ident)+".dat")):
	for line in open(path0+"single_particle/param_"+str(ident)+".dat","r"):
	        cont = line.split()
	        if(cont[0] == 'microstates'):
	                ms = int(cont[1])

	count =0 	
	W = np.zeros((bcount,ms,ms))
	Inv = np.zeros((bcount,ms,ms))
	for i in range(bmin,bmax):
		origin = str(ident + i*100)
		filename = path0+"single_particle/"+str(l)+"/reweight/T/Invariant_"+origin+"_"+aim+".dat"
		if (os.path.isfile(filename)):
			tran = np.loadtxt(filename)
			print(filename)
			Inv[i] = tran
			data = np.loadtxt(path0+"single_particle/count_"+ origin +".dat", dtype='i',usecols = [col])
			temp = (data[0:ms*ms]).reshape((ms,ms))
			W[i] = temp#(W[i] *count + temp ) / (count+ 1.)
			count +=1
		else:
			print("ALARM!! ALLAAAAAAARM!!")	


	I_opt = np.zeros((ms,ms))
	norm = np.zeros((ms,ms))
	
	for b in range(bmin,bmax):
		#if (np.count_nonzero(T[b]) > 0):
		I_opt = I_opt + W[b] * Inv[b]
		norm += W[b]


	for i in range(ms):
		for a in range(ms):
			if (norm[i,a] != 0):
				I_opt[i,a] /= norm[i,a]
	np.savetxt(path0+"single_particle/"+str(l)+"/reweight/T/Inv_opt_"+aim+".dat",I_opt)



	
