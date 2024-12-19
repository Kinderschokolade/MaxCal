import argparse
import numpy as np
import math
import getpass
import sys
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
from fct_markov import calc_av
from define_flow import def_flow_direction



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
wline = args.l
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



if (wline != None):
	wline= int( wline)

string = ' '
ms = 0
dT = 0
minpos = [0,0]
for line in open(path0+"Tgrad/param_"+str(ident)+".dat","r"):
	cont = line.split()
	if(cont[0] == 'dT'):
		dT = float(cont[1])
	if(cont[0] == 'microstates'):
		ms = int(cont[1])
	if(cont[0] == 'T'):
		Temperature1 = int(cont[1])
		Temperature2 = int(cont[2])
	if(cont[0] == 'xbas0'):
		minpos[0] = float(cont[1])
	if(cont[0] == 'xbas1'):
		minpos[1] = float(cont[1])
	if(cont[0] == 'gamma'):
		gamma = float(cont[1])
		

Temperature =  np.zeros(ms)
for i in range(ms):
	Temperature[i] = Temperature1 + (-Temperature1+ Temperature2) * i/(ms-1);

print(Temperature)



minpos[0] = math.floor( minpos[0] * ms) # for now only 1dim with 2 basins
minpos[1] = int(math.floor( minpos[1] * ms)) # for now only 1dim with 2 basins
print(minpos)

c_s = np.zeros(ms*ms)
p_s = np.zeros(ms)

F = open(path0+"Tgrad/count_"+ str(ident) +".dat","r") 
line = F.readline()
lagtime = line.split()
lagtime.pop(0)
lagtime  = [float(i) for i in lagtime]

cols = len(lagtime)

data = np.loadtxt(path0+"Tgrad/count_"+ str(ident) +".dat", dtype='i',usecols = (range(1,cols+1)))


###########
potential = np.loadtxt(path0+"Tgrad/potential_ms_"+str(ident)+".dat")
###############



c = data[0:ms*ms]
p = data[ms*ms:ms*ms+ms]   

L = np.zeros((ms,ms))
L = np.zeros((ms,ms))

lt = np.zeros((ms+1,cols))

EV_arr = np.zeros((ms+1,cols))
EV_comp = np.zeros((ms,2))
Evec_comp = np.zeros((ms,2))
Evecr_comp = np.zeros((ms,2))
Evec_arr = np.zeros((cols,ms,Enum))
Evecr_arr = np.zeros((cols,ms,Enum))
	

#print(c.shape[1])

T_check = np.zeros((ms,ms))

lvl = len(minpos)

mfpt_ar = np.zeros((lvl,lvl))
mfpt_arx = np.zeros((ms,ms))
mfpt_meas = np.zeros((ms,ms))
mfpt_err = np.zeros((ms,ms))
mfpt_steps = np.zeros((ms,ms))

F = def_flow_direction(ms)

blocks = 10
tra_len = 1000000

#r = flow_direction(cpl[:,wline].reshape((ms,ms)), cmi[:,wline].reshape((ms,ms)), ms)


tau = wline -1 ##FIND A BETTER WAY!!!

tao_MFPT = np.zeros((c.shape[1],10))
tao_J = np.zeros((c.shape[1],2))

for i in range(c.shape[1]):
	T = np.transpose(c[:,i].reshape((ms,ms)))
	q = p[:,i]
	
	for k in range(ms):
		for l in range(ms):
			if (q[k] > 0.):
				L[k,l] = T[k,l] / q[k]	
			else:
				L[k,l] = 0.
				L[k,k] = 1. 
	L[29,29] = 0.999;
	L[29,28] = 0.001
	L[28,29] = 0.001
	L[28,28] -= 0.001
	L[0,0] = 0.999
	L[0,1] =  0.001
	L[1,0] =  0.001
	L[1,1] -=  0.001
	Ev,Evec = np.linalg.eig(np.transpose(L))
	idx = Ev.argsort()[::-1] 
	Ev = Ev[idx]  #order eigenvectors by size3
	Evec = np.transpose(Evec)
	Evec = Evec[idx,:] 
	p_stat = np.real(Evec[0,:] / sum(Evec[0,:]))
	tao = lagtime[i] *dT
	lt[0,i] = tao	
	EV_arr[0,i] = tao # first line of Evals is tao (for gnuplot)
	tao_MFPT[i][0] = tao

#	if(np.linalg.cond(L) < 1./sys.float_info.epsilon):
#		for k in range(len(minpos)):
##			target = list(range(minpos[k] - rancut , minpos[k] + rancut +1))
#			for t in range(len(minpos)):
###				origin = list(range(minpos[t] - rancut , minpos[t] + rancut +1))
##				if(origin[-1] >= ms or target[-1] >= ms):
#					print ("origin/target for mfpt out of range!") 
#				else:
#					tao_MFPT[i][1+t*3+k] = mfpt(L,target = target, origin = origin)
#	else:
#		print(i, "matrix is singular")

	
	tao_J[i,0] = tao
	tao_J[i,1] = calc_av(F, L, p_stat)

	if (i == tau):	
		T_check = T
		L_check = L
		np.savetxt(path0+'Tgrad/'+str(int(tau))+'/EV/ps_'+str(ident)+'.dat',q/sum(q))
		p_check = q/sum(q)
		EV_comp[:,0] =  Ev.real
		EV_comp[:,1] =  Ev.imag

		# mean first passage time
		for k in range(len(minpos)):
			target = list(range(minpos[k] - rancut , minpos[k] + rancut +1))
			for t in range(len(minpos)):
				origin = list(range(minpos[t] - rancut , minpos[t] + rancut +1))
				if(origin[-1] >= ms or target[-1] >= ms):
					print ("origin/target for mfpt out of range!") 
				else:
					tau = wline * 1. - 1.  # read from first line count...
					mfpt_ar[t][k] = mfpt(L,target = target, origin = origin)

		#MSM_ = MSM(L)
		
		length = 1000
		fpt_ana = np.zeros((lvl*lvl,length))

		# First passage Time
		for k in range(lvl):
			for t in range(lvl):
				if (k!=t):
					start = list(range(minpos[k]-rancut, minpos[k]+rancut+1))
					end = list(range(minpos[t]-rancut, minpos[t]+rancut+1))
					fpt_ana[k*lvl+t,:] = first_passage(start,end,L,length)
					#print(start,end)		

		np.savetxt(path0+'Tgrad/'+str(int(tau))+'/MFPT/markovana_'+str(ident)+'.dat',np.transpose(fpt_ana))

		ratio = np.zeros((ms,ms))
		Sp = np.zeros((ms,ms))
		for k in range(ms):
			for l in range(ms):
				cutoff = 0.0001
				if (L[k,l]> cutoff and L[l,k] >cutoff):
					ratio[k,l] = np.log( L[k,l] / L[l,k])
			#	if (np.abs(k-l) <15):
#				force = np.sqrt(2*gamma*Temperature[k])- np.sqrt(2*gamma*Temperature[l])
#				Sp[k,l] = ( force / ms * (k-l)+ potential[k,1] - potential[l,1] ) / (Temperature[k])
				steps = abs(k-l)
				for a in range(steps):
					Sp[k,l] += 2*(potential[k,1] - potential[l,1])/ (Temperature[k] + Temperature[l])


			#	elif (l>k):
			#		force = np.sqrt(2*gamma*Temperature1)- np.sqrt(2*gamma*Temperature2)
			#		Sp[k,l] = ( force / ms * (ms-l+k)+ potential[k,1] - potential[l,1] ) / (Temperature1 + Temperature2)/2.
			#	else:
			#		force = np.sqrt(2*gamma*Temperature1)- np.sqrt(2*gamma*Temperature2)
			#		Sp[k,l] = ( -force / ms * (ms-k+l)+ potential[k,1] - potential[l,1] ) / (Temperature1 + Temperature2)/2.


		np.savetxt("/data/isilon/bause/Tgrad/"+str(int(tau))+"/T/ratio_"+str(ident)+".dat",ratio)
		np.savetxt("/data/isilon/bause/Tgrad/"+str(int(tau))+"/T/ratio_th_"+str(ident)+".dat",Sp)



		# test FPT by comparing to trajectory prododuced by MSM
#		tra = MSM_.simulate(tra_len)
		#stuff = mfpt_trajectory_ms_cross(tra, blocks, minpos, cut, lvl, rancut, ms,r)

#		stuff = mfpt_trajectory_ms(tra, blocks, minpos, cut, lvl, rancut, ms)
#		MFPT_hist = np.zeros((9,1000))
#		MFPT_meas = stuff[0]
#		MFPT_err = stuff[1]
#		MFPT_steps = stuff[2]	


#		np.savetxt('/data/isilon/bause/single_particle/MFPT/markov_'+str(int(tau))+"_"+str(ident)+'.dat',(mfpt_arx))
#		np.savetxt('/data/isilon/bause/single_particle/MFPT/measma_'+str(int(tau))+"_"+str(ident)+'.dat',(MFPT_meas))
#		np.savetxt('/data/isilon/bause/single_particle/MFPT/measmaerr_'+str(int(tau))+"_"+str(ident)+'.dat',(MFPT_err))
		


	

#		stuff = mfpt_trajectory_area(tra, blocks, minpos, cut, lvl, rancut, ms)

		MFPT_meas_reg = np.zeros((lvl,lvl))
		MFPT_err_reg = np.zeros((lvl,lvl))
		MFPT_steps_reg = np.zeros((lvl,lvl))
		MFPT_hist = np.zeros((lvl*lvl,1000))
#		MFPT_meas_reg = stuff[0]
#		MFPT_err_reg = stuff[1]
#		MFPT_steps_reg = stuff[2]	
#		MFPT_hist = stuff[3]	
#		for j in range(lvl):
#			for k in range(lvl):
#				if (j !=k ):
#					MFPT_hist[j*lvl+k,:] = MFPT_hist[j*lvl+k,:] /sum(MFPT_hist[j*lvl+k,:])
		print(tau, ident)
		np.savetxt(path0+'Tgrad/'+str(int(tau))+'/MFPT/markov_area_'+str(ident)+'.dat',(mfpt_ar))
#		np.savetxt('/data/isilon/bause/single_particle/MFPT/measma_area_'+str(int(tau))+"_"+str(ident)+'.dat',(MFPT_meas_reg))
#		np.savetxt('/data/isilon/bause/single_particle/MFPT/measmaerr_area_'+str(int(tau))+"_"+str(ident)+'.dat',(MFPT_err_reg))
#		np.savetxt('/data/isilon/bause/single_particle/MFPT/markovhist_'+str(int(tau))+"_"+str(ident)+'.dat',(np.transpose(MFPT_hist)))

		np.savetxt(path0+'Tgrad/'+str(int(tau))+'/T/T_'+str(ident)+'.dat',L)
		np.savetxt(path0+'Tgrad/'+str(int(tau))+'/EV/EV_com_'+str(ident)+'.dat',EV_comp)
		for k in range(Enum):
			Evec_comp[:,0] = Evec[k,:].real
			Evec_comp[:,1] = Evec[k,:].imag
			for l in range(ms):
				if (Evec[0,l] > 0):
					Evecr_comp[l,0] = (Evec[k,l]/ Evec[0,l]).real
					Evecr_comp[l,1] = (Evec[k,l]/ Evec[0,l]).imag
#				else:
#					Evecr_comp[l,0] = 0.
#					Evecr_comp[l,1] = 0.
			np.savetxt(path0+'Tgrad/'+str(int(tau))+'/EV/Evec'+str(k+1)+'_com_'+str(ident)+'.dat',Evec_comp)
			np.savetxt(path0+'Tgrad/'+str(int(tau))+'/EV/Evecr'+str(k+1)+'_com_'+str(ident)+'.dat',Evecr_comp)


	for j in range(1,ms):
		if (Ev[j].real > 0 and not math.isclose(Ev[j].real,1.0)): #How to deal with complex part?
			lt[j,i] = - tao / math.log(Ev[j].real)
		else:
			lt[j,i] = 'NaN'
		EV_arr[j+1,i] = Ev[j].real #fill array with Ev's in order
	
	for k in range(Enum):
		Evec_arr[i,:,k] = Evec[k+1].real # Evec[0] cotains invariant distr.
		for l in range(ms):
			if (Evec[0,l].real > 0):
				Evecr_arr[i,l,k] = Evec[k+1,l].real / Evec[0,l].real #right Evec
			else:
				Evecr_arr[i,l,k] = 0.		


np.savetxt(path0+"Tgrad/lagtime/J_"+ident+".dat",tao_J)	
np.savetxt(path0+"Tgrad/lagtime/MFPT_"+ident+".dat",tao_MFPT)	

if db_flag: # do det bal test if wanted
	db_mat = np.zeros((ms,ms))
	av_violation = 0.
	T_check = T_check.astype(float)
	
	T_check = T_check /  np.sum(T_check[:,:])	# normalise

	for i in range(ms):
		for j in range(ms):
			db_mat[i,j] = abs(T_check[i,j] -T_check[j,i])
			av_violation += db_mat[i,j]
	av_violation /= (ms*ms - ms)/2

	print("average violation of db: " ,av_violation)

	np.savetxt(path0+'Tgrad/'+str(int(tau))+'/detbal/err_'+str(ident)+'.dat',db_mat)
	np.savetxt(path0+'Tgrad/'+str(int(tau))+'/detbal/flux_'+str(ident)+'.dat',T_check)
	

if global_flag: # do global bal test if wanted
	err = np.zeros(ms)
	glob = np.zeros(ms)

	for i in range(ms):
		err[i] = abs( sum(p_check[j] * L_check[j,i] for j in range(ms) ) - p_check[i])
		glob[i] = abs( sum(p_check[j] * L_check[j,i] for j in range(ms) ))

	av_violation = sum(err) / ms
		
	p_check = np.transpose([p_check])
	glob = np.transpose([glob])

	out = np.concatenate((glob, p_check),axis = 1)

	np.savetxt(path0+'Tgrad/'+str(int(tau))+'/detbal/glob_'+str(ident)+'.dat',out)
	print("average violation of global balance: " ,av_violation)


np.savetxt(path0+'Tgrad/lagtime/'+str(ident)+'.dat',np.transpose(lt))
np.savetxt(path0+'Tgrad/EV/EV_'+str(ident)+'.dat',np.transpose(EV_arr))



for k in range(Enum):
	np.savetxt(path0+'Tgrad/EV/Evec'+str(k+1)+'_'+str(ident)+'.dat',np.transpose(Evec_arr[:,:,k]))

for k in range(Enum):
	np.savetxt(path0+'Tgrad/EV/Evecr'+str(k+1)+'_'+str(ident)+'.dat',np.transpose(Evecr_arr[:,:,k]))












