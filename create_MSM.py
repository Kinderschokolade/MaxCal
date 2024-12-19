import argparse
import numpy as np
import math
import getpass
from define_flow import read_minpos
from define_flow import read_cut
from msmtools.analysis import mfpt
#from pyemma.msm import MSM
#from analyse_trajectory import mfpt_trajectory_ms_cross
#from analyse_trajectory import mfpt_trajectory_area
#from analyse_trajectory import mfpt_trajectory_ms
#from define_cross import define_cross
#from define_flow import define_flow_direction
from fpt_ana import first_passage
from fct_markov import calc_Sprod
from fct_markov import calc_av
from define_flow import define_flow_direction
from scipy.optimize import fsolve
from fct_markov import moment_distribution
import pyemma.msm as msm
import msmtools 

def equ(c,ms,W):
	tup = [ (sum([W[j,k] *c[k]*c[j] for k in range(ms)] )-1.0) for j in range(ms)]
	tup = tuple(tup)
	return tup


def construct_MSM(T, ms, Sprod):
	W = [ [  np.sqrt(T[i,j] *T[j,i]) * math.exp(+ Sprod[i,j]/2. ) for j in range(ms)] for i in range(ms)]
	W = np.asarray(np.real(W))

	c = np.ones(ms)

	copt = fsolve(equ,c,args=(ms,W,),xtol=1e-8) # should converge for any starting point

	k = [[copt[j] *copt[i] *W[i,j] for j in range(ms)] for i in range(ms)]
	k= np.asarray(k)

	for i in range(ms): # minimal correction from numerics, eta = 1.
		k[i,:] = k[i,:] / sum(k[i,:]) 
	return np.real(k)


def symm_MSM(T, ms):
#	T = MSM(T,reversible=True)	
	T = msmtools.estimation.transition_matrix(T, reversible=True)
	return T



parser = argparse.ArgumentParser()
parser.add_argument('-o', help='ident')
parser.add_argument('-l', help='write T of specific line')
parser.add_argument('-n', help='number of (slowest) eigenvectors', default=3)
parser.add_argument('-db', help='set 1 for detbal check', default=0)
parser.add_argument('-gb', help='set 1 for global balance check', default=0)
parser.add_argument('-rc', help='rancut', default=1)
args = parser.parse_args();
ident = args.o
Enum = int(args.n)
tau = int(args.l)
db_flag = int(args.db)
global_flag = int(args.gb)
rancut = int(args.rc)

user = getpass.getuser()

if (user=="mbause"):
	path0 = "/u/mbause/data/"
	path1 = "/u/mbause/data/"
else:
	path0 = "/data/isilon/bause/"
	path1 = "/data/pckr194/bause/"



string = ' '
ms = 0
dT = 0
for line in open(path0+"single_particle/param_"+str(ident)+".dat","r"):
	cont = line.split()
	if(cont[0] == 'dT'):
		dT = float(cont[1])
	if(cont[0] == 'microstates'):
		ms = int(cont[1])
	if(cont[0] == 'T0'): # old version
		Temperature = float(cont[1])
	if(cont[0] == 'T'): # new version
		Temperature = float(cont[1])
	if(cont[0] == 'extf'):
		force = float(cont[1])


c_s = np.zeros(ms*ms)
p_s = np.zeros(ms)

F = open(path0+"single_particle/count_"+ str(ident) +".dat","r") 
line = F.readline()
lagtime = line.split()
lagtime.pop(0)
lagtime  = [int(i) for i in lagtime]

cols = len(lagtime)

data = np.loadtxt(path0+"single_particle/count_"+ str(ident) +".dat", dtype='i',usecols = (range(1,cols+1)))

###########
potential = np.loadtxt(path0+"single_particle/potential_ms_"+str(ident)+".dat")
###############

ms, minpos = read_minpos(path0,ident)

c = data[0:ms*ms]
p = data[ms*ms:ms*ms+ms]   

Sprod = np.loadtxt(path0+"single_particle/Sprod_meas_"+str(ident)+".dat")

L = np.zeros((ms,ms))

EV_arr = np.zeros((ms+1,cols))
EV_comp = np.zeros((ms,2))
Evec_comp = np.zeros((ms,2))
Evecr_comp = np.zeros((ms,2))
Evec_arr = np.zeros((cols,ms,Enum))
Evecr_arr = np.zeros((cols,ms,Enum))
	
lvl = len(minpos)
lt = np.zeros((ms+1,cols))

mfpt_ar = np.zeros((lvl,lvl))
mfpt_meas = np.zeros((ms,ms))
mfpt_err = np.zeros((ms,ms))
mfpt_steps = np.zeros((ms,ms))

tao_MFPT = np.zeros((c.shape[1],10))
tao_J = np.zeros((c.shape[1],2))
tao_Sprod = np.zeros((c.shape[1],2))

steps = 100;
p_ck = np.zeros((lvl+1,steps))
p_tr = np.zeros((2*lvl+1,c.shape[1]))

offequ = np.zeros((ms,ms))
#escape time trajectory based (ck-test)

for i in range(cols):
	T = np.transpose(c[:,i].reshape((ms,ms)))
	q = p[:,i]
	for k in range(ms):
		for l in range(ms):
			if(q[k] > 0.):
				L[k,l] = T[k,l] / sum(T[k,:]) #/ q[k]	
			else:
				L[k,l] = 0.
			#	L[k,k] = 1. 



	Ev,Evec = np.linalg.eig(np.transpose(L))
	idx = Ev.argsort()[::-1] 
	Ev = Ev[idx]  #order eigenvectors by size
	Evec = np.transpose(Evec)
	Evec = Evec[idx,:] 
	p_stat = np.real(Evec[0,:] / sum(Evec[0,:]))
	tao = lagtime[i]
	lt[0,i] = tao	
	EV_arr[0,i] = tao # first line of Evals is tao (for gnuplot)
	tao_MFPT[i][0] = tao
	
	p_tr[0,i] = tao;

	tao_J[i,0] = tao
	tao_Sprod[i,0] = tao
	tao_J[i,1] = 0. #calc_av(F, L, p_stat)  # error because p_stat is of antoher size
	tao_Sprod[i,1] = calc_Sprod(L)

	if (lagtime[i] == tau):	
		L_check = L
		L_symm= symm_MSM(T,ms)
		L_con=construct_MSM(L, ms, Sprod)
		#np.savetxt(path0+'single_particle/'+str(int(tau))+'/EV/ps_'+str(ident)+'.dat',q/sum(q))
		p_check = q/sum(q)
			
		EV_comp[:,0] =  Ev.real # size problem
		EV_comp[:,1] =  Ev.imag

		Tc = np.transpose(c[:,i].reshape((ms,ms)))
		T = np.zeros((ms,ms))
		# mean first passage time
		for k in range(len(minpos)):
			target = list(range(minpos[k] - rancut , minpos[k] + rancut +1))
			for t in range(len(minpos)):
				origin = list(range(minpos[t] - rancut , minpos[t] + rancut +1))
				if(origin[-1] >= ms or target[-1] >= ms):
					print ("origin/target for mfpt out of range!") 
#				else:
#					tau = wline * 1. - 1.  # read from first line count...
#					mfpt_ar[t][k] = mfpt(L,target = target, origin = origin)

		ratio = np.zeros((ms,ms))
		ratio_con = np.zeros((ms,ms))
		ratio_symm = np.zeros((ms,ms))
		for k in range(ms):
			for l in range(ms):
				if (L[k,l] > 0. and L[l,k] > 0. ):
					ratio[k,l] = L[k,l] / L[l,k] 
					ratio_con[k,l] = L_con[k,l] / L_con[l,k] 
					ratio_symm[k,l] = L_symm[k,l] / L_symm[l,k] 

		length = 10000
		fpt_ana = np.zeros((lvl*lvl,length))
		fpt_ana_con = np.zeros((lvl*lvl,length))
		fpt_ana_symm = np.zeros((lvl*lvl,length))
		mom = np.zeros((lvl*lvl,3))
		mom_con = np.zeros((lvl*lvl,3))
		mom_symm = np.zeros((lvl*lvl,3))

		# First passage Time
		for k in range(lvl):
			for t in range(lvl):
				if (k!=t):
					start = list(range(minpos[k]-rancut, minpos[k]+rancut+1))
					end = list(range(minpos[t]-rancut, minpos[t]+rancut+1))
					fpt_ana[k*lvl+t,:] = first_passage(start,end,L,length)
					fpt_ana_con[k*lvl+t,:] = first_passage(start,end,L_con,length)
					fpt_ana_symm[k*lvl+t,:] = first_passage(start,end,L_symm,length)
					for m in range(1,4):
						mom_con[k*lvl+t,m-1] = moment_distribution(m,fpt_ana_con[k*lvl+t])
						mom[k*lvl+t,m-1] = moment_distribution(m,fpt_ana[k*lvl+t])
						mom_symm[k*lvl+t,m-1] = moment_distribution(m,fpt_ana_symm[k*lvl+t])

					mom[k*lvl+t,1] = np.sqrt( mom[k*lvl+t,1] - mom[k*lvl+t,0]* mom[k*lvl+t,0] ) #standard deviation
					mom[k*lvl+t,2] = (mom[k*lvl+t,2] - 3.*mom[k*lvl+t,0] * mom[k*lvl+t,1] *mom[k*lvl+t,1]  - mom[k*lvl+t,0]* mom[k*lvl+t,0]* mom[k*lvl+t,0])/ (mom[k*lvl+t,1] * mom[k*lvl+t,1] * mom[k*lvl+t,1] ) #skewness
					mom_con[k*lvl+t,1] = np.sqrt( mom_con[k*lvl+t,1] - mom_con[k*lvl+t,0]* mom_con[k*lvl+t,0] ) #standard deviation
					mom_con[k*lvl+t,2] = (mom_con[k*lvl+t,2] - 3.*mom_con[k*lvl+t,0] * mom_con[k*lvl+t,1] *mom_con[k*lvl+t,1]  - mom_con[k*lvl+t,0]* mom_con[k*lvl+t,0]* mom_con[k*lvl+t,0])/ (mom_con[k*lvl+t,1] * mom_con[k*lvl+t,1] * mom_con[k*lvl+t,1] ) #skewness
					mom_symm[k*lvl+t,1] = np.sqrt( mom_symm[k*lvl+t,1] - mom_symm[k*lvl+t,0]* mom_symm[k*lvl+t,0] ) #standard deviation
					mom_symm[k*lvl+t,2] = (mom_symm[k*lvl+t,2] - 3.*mom_symm[k*lvl+t,0] * mom_symm[k*lvl+t,1] *mom_symm[k*lvl+t,1]  - mom_symm[k*lvl+t,0]* mom_symm[k*lvl+t,0]* mom_symm[k*lvl+t,0])/ (mom_symm[k*lvl+t,1] * mom_symm[k*lvl+t,1] * mom_symm[k*lvl+t,1] ) #skewness
	
		print(mom[3,0],mom_symm[3,0], mom_con[3,0])

		for k in range(1,4):
			np.savetxt("/data/isilon/bause/single_particle/"+str(int(tau))+"/MFPT/moment_con"+str(int(k))+"_"+str(ident)+".dat",mom_con[:,k-1].reshape((lvl,lvl)))
			np.savetxt("/data/isilon/bause/single_particle/"+str(int(tau))+"/MFPT/moment_symm"+str(int(k))+"_"+str(ident)+".dat",mom_symm[:,k-1].reshape((lvl,lvl)))
			np.savetxt("/data/isilon/bause/single_particle/"+str(int(tau))+"/MFPT/moment"+str(int(k))+"_"+str(ident)+".dat",mom[:,k-1].reshape((lvl,lvl)))
	
		np.savetxt(path0+'single_particle/'+str(int(tau))+'/MFPT/markovana_con_'+str(ident)+'.dat',np.transpose(fpt_ana_con))
		np.savetxt(path0+'single_particle/'+str(int(tau))+'/MFPT/markovana_symm_'+str(ident)+'.dat',np.transpose(fpt_ana_symm))
		np.savetxt(path0+'single_particle/'+str(int(tau))+'/MFPT/markovana_'+str(ident)+'.dat',np.transpose(fpt_ana))
		np.savetxt(path0+'single_particle/'+str(int(tau))+'/T/ratio_con_'+str(ident)+'.dat',ratio_con)
		np.savetxt(path0+'single_particle/'+str(int(tau))+'/T/ratio_symm_'+str(ident)+'.dat',ratio_symm)
		np.savetxt(path0+'single_particle/'+str(int(tau))+'/T/ratio_'+str(ident)+'.dat',ratio)

		ratio = np.zeros((ms,ms))
		act = np.zeros((ms,ms))
		ratio_err =0.
		#for k in range(ms):
		#	for l in range(ms):
#
#				cutoff = 0.00
#				act[k,l] = L_real[k,l] * L_real[l,k] 
#				if (L_real[k,l]> cutoff and L_real[l,k] >cutoff):
#					ratio[k,l] =  L_real[k,l] / L_real[l,k]
#					ratio_err += L_real[k,l]* np.abs(Sp[k,l]-ratio[k,l])/ms
#				else:
#					ratio[k,l]= np.nan 
						#print(Sp[k,l],ratio[k,l], np.log(ratio[k,l]/Sp[k,l]))
		
#		ratio_err /= maxi
#		print("ratio error",ratio_err)


#		np.savetxt("/data/isilon/bause/single_particle/"+str(int(tau))+"/T/ratio_"+str(ident)+".dat",ratio)
#		np.savetxt("/data/isilon/bause/single_particle/"+str(int(tau))+"/T/act_"+str(ident)+".dat",act)
#		np.savetxt("/data/isilon/bause/single_particle/"+str(int(tau))+"/T/ratio_th_"+str(ident)+".dat",np.exp(Sp))
#		np.savetxt("/usr/data/bause/single_particle/"+str(int(tau))+"/T/Sprod_meas_"+str(ident)+".dat",Sprod)
		np.savetxt(path0+'single_particle/'+str(int(tau))+'/EV/EV_com_'+str(ident)+'.dat',EV_comp)

		for k in range(Enum):
			Evec_comp[:,0] = Evec[k,:].real
			Evec_comp[:,1] = Evec[k,:].imag
			for l in range(ms):
				if (Evec[0,l] > 0):
					Evecr_comp[l,0] = (Evec[k,l]/ Evec[0,l]).real
					Evecr_comp[l,1] = (Evec[k,l]/ Evec[0,l]).imag
				else:
					Evecr_comp[l,0] = 0.
					Evecr_comp[l,1] = 0.
			np.savetxt(path0+'single_particle/'+str(int(tau))+'/EV/Evec'+str(k+1)+'_com_'+str(ident)+'.dat',Evec_comp)
			np.savetxt(path0+'single_particle/'+str(int(tau))+'/EV/Evecr'+str(k+1)+'_com_'+str(ident)+'.dat',Evecr_comp)

	msred = L.shape[0]
	for j in range(1,msred):
		if (Ev[j].real > 0 and not math.isclose(Ev[j].real,1.0)): #How to deal with complex part?
			lt[j,i] = - tao / math.log(Ev[j].real)
		else:
			lt[j,i] = 'NaN'
		EV_arr[j+1,i] = Ev[j].real #fill array with Ev's in order

np.savetxt(path0+"single_particle/lagtime/J_"+ident+".dat",tao_J)	
np.savetxt(path0+"single_particle/lagtime/ck_"+ident+".dat",np.transpose(p_tr))	
np.savetxt(path0+"single_particle/lagtime/Sprod_"+ident+".dat",tao_Sprod)	
np.savetxt(path0+"single_particle/lagtime/MFPT_"+ident+".dat",tao_MFPT)	


np.savetxt(path0+'single_particle/lagtime/'+str(ident)+'.dat',np.transpose(lt))
np.savetxt(path0+'single_particle/EV/EV_'+str(ident)+'.dat',np.transpose(EV_arr))



for k in range(Enum):
	np.savetxt(path0+'single_particle/EV/Evec'+str(k+1)+'_'+str(ident)+'.dat',np.transpose(Evec_arr[:,:,k]))

for k in range(Enum):
	np.savetxt(path0+'single_particle/EV/Evecr'+str(k+1)+'_'+str(ident)+'.dat',np.transpose(Evecr_arr[:,:,k]))












