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


def reweight_min(x,T,F,E,ms,J_comp,E_comp):
	err,k,p=reweight(x,T,F,E,ms,J_comp,E_comp)
	return err

def reweight(x,T,F,E, ms,J_comp, E_comp):
	#W = [ [ np.power( T[i,j] , 2./3.) * np.power(T[j,i], 1./3.) * math.exp(2./3.*gamma*r[i,j] + 1./3. *gamma*r[j,i]) for j in range(ms)] for i in range(ms)]
	W = [ [  T[i,j] * math.exp(x[0]*F[i,j]+x[1]*E[i,j]) for j in range(ms)] for i in range(ms)]
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

	Jerr = np.abs(J_comp - calc_av(F,k,p))
	Eerr = np.abs(E_comp - calc_av(E,k,p))
	err = Jerr+Eerr

	return err, np.real(k),np.real(p)



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

Ev,Evec =  analyse_MSM(T)
p_init = np.real(Evec[0] / sum(Evec[0]) )


Ev,Evec =  analyse_MSM(T_comp)
p_comp = np.real(Evec[0] / sum(Evec[0]) )
temp = np.loadtxt("/data/isilon/bause/single_particle/potential_ms_"+str(ident)+".dat")
Esur = np.zeros((ms,ms))

for i in range(ms):
	for j in range(ms):
		Esur[i,j] = (temp[i,1] - temp[j,1]) / 2.

sp= np.zeros((ms,ms))
ratio_t = np.zeros((ms,ms))
cutoff = 0.0000001
for k in range(ms):
	for l in range(ms):
		if (np.abs(k-l) <15):
				sp[k,l] = ( force / ms * (k-l)+ potential[k,1] - potential[l,1] ) / temperature
		elif (l>k):
				sp[k,l] = ( force / ms * (ms-l+k)+ potential[k,1] - potential[l,1] ) / temperature
		else:
				sp[k,l] = ( -force / ms * (ms-k+l)+ potential[k,1] - potential[l,1] ) / temperature

Esur = sp

lvl = len(minpos)
Ng = 800
gamma = np.zeros(Ng)
MFPT_ar = np.zeros((Ng, lvl *lvl))
ratio = np.zeros((Ng,3))

J_ar = np.zeros(Ng)

MFPT_comp = mfpt_markov(T_comp , minpos, rancut, ms,lvl).reshape(lvl*lvl)
ratio_comp =  sum( [MFPT_comp[1]/MFPT_comp[3] , MFPT_comp[6]/MFPT_comp[2] , MFPT_comp[5]/MFPT_comp[7]] ) / 3.
J_comp =  calc_av(r, T_comp, p_comp)
E_comp =  calc_av(Esur, T_comp, p_comp)
E_init =  calc_av(Esur, T, p_init)


Sprod_comp = calc_Sprod(T_comp)

#ratio_cur = 10
J_cur = 0
J_err = 100

Jdiff = 1.
diff_min = 0.00001

out_var = fmin(reweight_min,x0=np.array([0.,0.]),ftol = diff_min , disp=False, args=(T ,r , Esur, ms , J_comp, E_comp,))
gamma = out_var[0]
delta = out_var[1]

J_err, T_cur, p_cur = reweight(out_var,T,r,Esur,ms,J_comp,E_comp)
J_cur = calc_av(r,T_cur,p_cur)

Sprod_rew = calc_Sprod(T_cur)
MFPT_rew = mfpt_markov(T_cur , minpos, rancut, ms,lvl).reshape(lvl*lvl)
ratio_rew =  np.sum([MFPT_rew[1]/MFPT_rew[3] , MFPT_rew[6]/MFPT_rew[2] , MFPT_rew[5]/MFPT_rew[7]]) / 3.
E_cur = calc_av(Esur, T_cur, p_cur)


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
	np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/MFPT/moment"+str(int(k))+"_"+str(identin)+"_"+str(ident)+".dat",mom_cur[:,k-1].reshape((lvl,lvl)))
	np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/MFPT/moment"+str(int(k))+"_"+str(ident)+".dat",mom_comp[:,k-1].reshape((lvl,lvl)))


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
for k in range(ms):
	for l in range(ms):
		if (T_cur[k,l] > cutoff and T_cur[l,k] > cutoff):
			ratio_t[k,l] = np.log( T_cur[k,l] / T_cur[l,k] )


np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/T/ratio_"+str(identin)+"_"+str(ident)+".dat",ratio_t)
np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/T/ratio_th_"+str(identin)+"_"+str(ident)+".dat",sp)
np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/T/ratio_diff_"+str(identin)+"_"+str(ident)+".dat",np.abs(sp-ratio_t))
	
p_cur = np.transpose([p_cur])
p_comp = np.transpose([p_comp])
out = np.concatenate((p_cur, p_comp),axis = 1)
np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/T/p_"+str(ident)+".dat", out)




print(identin, ident,gamma,delta)
print(J_cur,J_comp)
print(E_cur, E_comp, E_init)
print(Sprod_comp, Sprod_rew)


MFPT_cumhist = hist_to_cumhist(MFPT_hist)
#print(identin, ident,gamma_cur, tot_error)
np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/MFPT/hist_"+str(identin)+"_"+str(ident)+".dat", np.transpose(MFPT_hist))
np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/MFPT/cumhist_"+str(identin)+"_"+str(ident)+".dat", np.transpose(MFPT_cumhist))
np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/T/T_"+str(identin)+"_"+str(ident)+".dat", T_cur)

	
#gamma = np.transpose([gamma])
#J_ar = np.transpose([J_ar])
#out = np.concatenate((gamma, J_ar),axis = 1)
#np.savetxt("/data/isilon/bause/single_particle/"+str(lag)+"/reweight/J/J_"+str(identin)+".dat", out)




