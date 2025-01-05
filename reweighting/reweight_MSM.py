import math
import argparse
import numpy as np
import os
from define_flow import define_flow_direction
from define_flow import read_cut
from define_flow import pos_pot

def calc_err_relaxtime_r(T, r, Ev_comp, gamma):
	N = len(Ev_comp)

	W = [ [ T[i,j] * math.exp( gamma * r[i,j] ) for j in range(ms)] for i in range(ms)]
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
	p = Evecl[0,:] * Evecr[0,:]
	p = np.asarray(p)
	
	p = p / sum(p)

	# may cause divide by 0 error	
	k = [[ Evecr[0,j] / (Evecr[0,i] * Ev[0]) *W[i][j] for j in range(ms)] for i in range(ms)]
	
	#find relaxation times of k	
	Ev,Evec = np.linalg.eig(np.transpose(k)) # left eigenvector
	idx = Ev.argsort()[::-1] 
	Ev = Ev[idx]  #order eigenvectors by size
	Evec = np.transpose(Evec)
	Evec = Evec[idx,:] 
	#print(p/sum(p))
	#print(gamma,Ev[:N])
	rel = np.zeros(N)
	rel_comp = np.zeros(N)
	err = 0.
	for i in range(1,N):
		rel[i] = 1./math.log(np.real(Ev[i]))
		rel_comp[i] = 1./math.log(np.real(Ev_comp[i]))
		err += (rel[i] - rel_comp[i])*(rel[i] - rel_comp[i]) 
	return err



def calc_err_relaxtime(Tmi, Tpl, Ev_comp, gamma):
	N = len(Ev_comp)
	W = [ [ Tmi[i,j] * math.exp( gamma ) + Tpl[i,j] * math.exp( -gamma) for j in range(ms)] for i in range(ms)]
	Ev,Evecl = np.linalg.eig(np.transpose(W))
	Ev,Evecr = np.linalg.eig(W)
	p = Evecl[:,0] * Evecr[:,0]
	p = np.asarray(p)
	p = p / sum(p)
	k = [[ Evecr[j,0] / (Evecr[i,0] * Ev[0]) *W[i][j] for j in range(ms)] for i in range(ms)]
	
	#find relaxation times of k	
	Ev,Evec = np.linalg.eig(np.transpose(k)) # left eigenvector
	idx = Ev.argsort()[::-1] 
	Ev = Ev[idx]  #order eigenvectors by size
	Evec = np.transpose(Evec)
	Evec = Evec[idx,:] 
	#print(p/sum(p))
	#print(gamma,Ev[:N])
	rel = np.zeros(N)
	rel_comp = np.zeros(N)
	err = 0.
	for i in range(1,N):
		rel[i] = 1./math.log(np.real(Ev[i]))
		rel_comp[i] = 1./math.log(np.real(Ev_comp[i]))
		err += (rel[i] - rel_comp[i])*(rel[i] - rel_comp[i]) 
	return err

def calc_err_flux(Tmi, Tpl, flux_comp, gamma):
	ms = Tmi.shape[0]
	W = [ [ Tmi[i,j] * math.exp( gamma ) + Tpl[i,j] * math.exp( -gamma) for j in range(ms)] for i in range(ms)]
	Ev,Evecl = np.linalg.eig(np.transpose(W))
	Ev,Evecr = np.linalg.eig(W)
	p = Evecl[:,0] * Evecr[:,0]
	p = np.asarray(p)
	p = p / sum(p)
	k = [[ Evecr[j,0] / (Evecr[i,0] * Ev[0]) *W[i][j] for j in range(ms)] for i in range(ms)]
	k = np.asarray(k)
	for i in range(ms):
		k[i] = k[i] / sum(k[i,:])
	
	fl = flux(p,k)
	err = (flux_comp - fl ) * (flux_comp - fl)
	print(gamma, flux_comp, fl)

	return err

def hist_reweight(identin, identout):
	N =2
	start = 6000
	stop = 6900
	fl = np.loadtxt("/data/isilon/bause/single_particle/reweight/weight_" + str(start)+"_"+str(stop) +".dat")
	i = int((int(identin) -start )/ 100)
	j = int((int(identout) -start) /100)
	#print(i,j, fl[i,j])
	gamma= fl[i,j]

	return gamma


def calc_err_Evec(Tmi,Tpl,Ev_comp,gamma):
	
	W = [ [ Tmi[i,j] * math.exp( gamma ) + Tpl[i,j] * math.exp( -gamma) for j in range(ms)] for i in range(ms)]
	Ev,Evecl = np.linalg.eig(np.transpose(W))
	Ev,Evecr = np.linalg.eig(W)
	p = Evecl[:,0] * Evecr[:,0]
	p = np.asarray(p)
	p = p / sum(p)
	k = [[ Evecr[j,0] / (Evecr[i,0] * Ev[0]) *W[i][j] for j in range(ms)] for i in range(ms)]
	
	#find relaxation times of k	
	Ev,Evec = np.linalg.eig(np.transpose(k)) # left eigenvector
	idx = Ev.argsort()[::-1] 
	Ev = Ev[idx]  #order eigenvectors by size
	Evec = np.transpose(Evec)
	Evec = Evec[idx,:] 
	
	err1 = sum((Evec[1,i] - Ev_comp[i,1])*(Evec[1,i] - Ev_comp[i,1]) for i in range(ms) )
	err2 = sum((Evec[1,i] + Ev_comp[i,1])*(Evec[1,i] + Ev_comp[i,1]) for i in range(ms) ) 
	# Evec could differ a minus sign
	err = min(err1,err2)
	#print(gamma, err)	
	return err


def calc_err_Evec_r(T,r,Ev_comp,gamma):
	
	W = [ [ T[i,j] * math.exp( gamma * r[i,j] ) for j in range(ms)] for i in range(ms)]
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
	p = Evecl[0,:] * Evecr[0,:]
	p = np.asarray(p)
	
	p = p / sum(p)

	# may cause divide by 0 error	
	k = [[ Evecr[0,j] / (Evecr[0,i] * Ev[0]) *W[i][j] for j in range(ms)] for i in range(ms)]
	

	#find relaxation times of k	
	Ev,Evec = np.linalg.eig(np.transpose(k)) # left eigenvector
	idx = Ev.argsort()[::-1] 
	Ev = Ev[idx]  #order eigenvectors by size
	Evec = np.transpose(Evec)
	Evec = Evec[idx,:] 
	
	err1 = sum((Evec[1,i] - Ev_comp[i,1])*(Evec[1,i] - Ev_comp[i,1]) for i in range(ms) )
	err2 = sum((Evec[1,i] + Ev_comp[i,1])*(Evec[1,i] + Ev_comp[i,1]) for i in range(ms) )
	err = min(err1,err2)

	#print(gamma, err)	
	return err



def flux(p, k):
	ms = p.shape[0]
	# this cannot work, because Tpl and Tmi is unknown after MSM reweighting!!
	return (sum(sum( p[i] * k[i,j] for j in range(ms)) for i in range(ms)))

parser = argparse.ArgumentParser()
parser.add_argument('-o', help='reweight to ident ')
parser.add_argument('-i', help='reweight from ident')
args = parser.parse_args();
ident = args.i
identout = args.o


N = 2 # number of eigenvalues to compare

# find Eigensystem for comparison
path = "/data/isilon/bause/single_particle/T/Tpl_"+ str(identout) +".dat"
Tpl = np.loadtxt(path)
path = "/data/isilon/bause/single_particle/T/Tmi_"+ str(identout) +".dat"
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

Ev_comp,Evec_comp = np.linalg.eig(np.transpose(Tpl+Tmi))

Ev_comp = Ev_comp[:N]

p = Evec_comp[:,0]/sum(Evec_comp[:,0])

gamma_min = 0.0


path = "/data/isilon/bause/single_particle/T/Tmi_"+ str(ident) +".dat"
Tmi = np.loadtxt(path)
path = "/data/isilon/bause/single_particle/T/Tpl_"+ str(ident) +".dat"
Tpl = np.loadtxt(path)
ms = np.shape(Tpl)[0]

for i in range(ms):
	denom = sum(Tmi[i,:]+Tpl[i,:])

	for j in range(ms):
		if (denom > 0.):
			Tmi[i,j] = Tmi[i,j] / denom
			Tpl[i,j] = Tpl[i,j] / denom
		else:
			Tmi[i,j] = 0.
			Tpl[i,j] = 0.




Tpl = np.asarray(Tpl,dtype = np.float32)
Tmi = np.asarray(Tmi,dtype = np.float32)

T = Tpl+Tmi

ms,cut = read_cut(ident)
r = define_flow_direction(ms,cut)

np.savetxt("test.dat", r)
err_min = calc_err_relaxtime_r(T,r, Ev_comp, gamma_min)
#err_min = calc_err_flux(Tmi,Tpl, flux_comp, gamma_min)
#err_min = calc_err_Evec_r(T,r,Evec_comp,gamma_min)


#gamma_hist = hist_reweight(ident,identout)

for i in range(-50,50):
	gamma =   i / 100. #random choice here
	#err = calc_err_flux(Tmi,Tpl, flux_comp, gamma)
	#err = calc_err_Evec_r(T,r,Evec_comp,gamma)
	err = calc_err_relaxtime_r(T,r,Ev_comp, gamma)
	if (err < err_min):
		err_min = err
		gamma_min =gamma
	print(gamma, err)

gamma= gamma_min

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
	
Ev,Evec = np.linalg.eig(np.transpose(k))
idx = Ev.argsort()[::-1] 
Ev = Ev[idx]  #order eigenvectors by size
Evec = np.transpose(Evec)
Evec = Evecl[idx,:] 

#print("weighted ",Ev[:N])
#print("reference ",Ev_comp)
#print(gamma_min, err_min)

print(ident, identout, gamma)


with open("/data/isilon/bause/single_particle/reweight/gamma.dat" , 'a') as file:
	file.write(ident+"\t"+identout +"\t"+ str(-1./math.log(np.real(Ev[1]))) + "\t" + str(-1./math.log(np.real(Ev_comp[1])))+ "\n")

np.savetxt("/data/isilon/bause/single_particle/reweight/"+str(ident)+"_"+ str(identout)+ ".dat",np.real(k))

out = np.transpose([np.asarray(Evec[:,1]), np.asarray(Evec_comp[:,1])])

np.savetxt("/data/isilon/bause/single_particle/reweight/Evec1_"+str(ident)+"_"+ str(identout)+ ".dat",np.real(out))


out = np.transpose([np.asarray(Evec[:,0]/sum(Evec[:,0])), np.asarray(Evec_comp[:,0]/sum(Evec_comp[:,0]))])
np.savetxt("/data/isilon/bause/single_particle/reweight/Evec0_"+str(ident)+"_"+ str(identout)+ ".dat",np.real(out))


np.savetxt("/data/isilon/bause/single_particle/reweight/lagtime_"+str(ident)+"_"+ str(identout)+ ".dat",-1./np.log(np.real(Ev[1:])))

