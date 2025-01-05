import argparse
import getpass
import math

import numpy as np
from define_flow import read_cut, read_minpos
from msmtools.analysis import mfpt

# from pyemma.msm import MSM
# from analyse_trajectory import mfpt_trajectory_ms_cross
# from analyse_trajectory import mfpt_trajectory_area
# from analyse_trajectory import mfpt_trajectory_ms
# from define_cross import define_cross
# from define_flow import define_flow_direction
from fpt_ana import first_passage

# uses -o count matrices (for different lagtime)
# time scale. Use to identify markovian region
# choose -l to fix lagtime and perform principal component analysis - Complex for off-equ!

parser = argparse.ArgumentParser()
parser.add_argument("-o", help="ident")
parser.add_argument("-l", help="write T of specific line")
parser.add_argument("-n", help="number of (slowest) eigenvectors", default=3)
parser.add_argument("-db", help="set 1 for detbal check", default=0)
parser.add_argument("-rc", help="rancut", default=1)
args = parser.parse_args()
ident = args.o
Enum = int(args.n) + 1
tau = int(args.l)
db_flag = int(args.db)
rancut = int(args.rc)

user = getpass.getuser()


# change path0 to the directory, where you want the data to be saved
# also use the same folder structure as we use in this script. (use python ( e.g. with 'import os') or batch to create all folders, it will take forever to do it by hand)
# if numpy.savetxt does not find the folder, it simply does not save the data.
if user == "mbause":
    path0 = "/u/mbause/data/"
    path1 = "/u/mbause/data/"
else:
    path0 = "/data/isilon/bause/"
    path1 = "/data/pckr194/bause/"


string = " "
MS = [0, 0]
dT = 0

# param file, created by the simulation - MS in x and y direction is predefined  here
for line in open(path0 + "single_particle/param_" + str(ident) + ".dat", "r"):
    cont = line.split()
    if cont[0] == "dT":
        dT = float(cont[1])
    if cont[0] == "microstates":
        MS[0] = int(cont[1])
        MS[1] = int(cont[3])

ms = MS[0] * MS[1]

c_s = np.zeros(ms * ms)
p_s = np.zeros(ms)

# count file contains count matrices for all lagtime (0th column - tao = 1 , 1st column - tao =1 , 2nd column - tao=2 ,...) , to  a max of 'cols'
F = open(path0 + "single_particle/count_" + str(ident) + ".dat", "r")
line = F.readline()
lagtime = line.split()
lagtime.pop(0)
lagtime = [int(i) for i in lagtime]
# alternatively you can create count and param yourself - by your choice of microstates or by applying an algorithm to create microstates.

cols = len(lagtime)
data = np.loadtxt(path0 + "single_particle/count_" + str(ident) + ".dat", dtype="i", usecols=(range(1, cols + 1)))

# minpos gives the position of the minimum in microstates
# cut is the dividing microstate between two minima
ms, minpos = read_minpos(path0, ident)
ms, cut = read_cut(path0, ident)
# This was used for 1D and not applicable anymore


ms = MS[0] * MS[1]
# the first ms*ms lines cotain countmatrix information
# the last ms lines contain density informaton
c = data[0 : ms * ms]
p = data[ms * ms : ms * ms + ms]

L = np.zeros((ms, ms))

lt = np.zeros((ms + 1, cols))

EV_arr = np.zeros((ms + 1, cols))
EV_comp = np.zeros((ms, 2))
Evec_comp = np.zeros((ms, 2))
Evecr_comp = np.zeros((ms, 2))
Evec_arr = np.zeros((cols, ms, Enum))
Evecr_arr = np.zeros((cols, ms, Enum))

T_check = np.zeros((ms, ms))

lvl = len(minpos)

mfpt_ar = np.zeros((lvl, lvl))
mfpt_arx = np.zeros((ms, ms))
mfpt_meas = np.zeros((ms, ms))
mfpt_err = np.zeros((ms, ms))
mfpt_steps = np.zeros((ms, ms))


pstat = np.zeros((MS[0], MS[1]))
for i in range(c.shape[1]):
    T = c[:, i].reshape((ms, ms))
    q = p[:, i]

    for k in range(ms):
        for l in range(ms):
            if q[k] > 0.0:
                L[k, l] = T[k, l] / q[k]  # transition matrix
            else:
                # L[k,l] = 0.
                # L[k,k] = 1.
                L = L.delete(L, (k), axis=0)
                L = L.delete(L, (k), axis=1)

    print(lagtime[i], L.shape)
    Ev, Evec = np.linalg.eig(np.transpose(L))  # Eigenvalue and Eigenvector
    idx = Ev.argsort()[::-1]
    Ev = Ev[idx]  # order eigenvalues by size of EV
    Evec = np.transpose(Evec)
    Evec = Evec[idx, :]
    tao = lagtime[i] * dT  # lagtime in seconds
    lt[0, i] = tao
    EV_arr[0, i] = tao  # first line of Evals is tao (for gnuplot)

    taus = str(int(tau))  # the string of lagtime
    if lagtime[i] == tau:  ## focus on tau
        print(tau)
        for k in range(MS[0]):
            for l in range(MS[1]):
                pstat[k, l] = q[l * MS[0] + k]
        T_check = T
        np.savetxt(
            path0 + "single_particle/" + taus + "/EV/ps_" + str(ident) + ".dat", pstat
        )  # measured steady state density
        np.savetxt(path0 + "single_particle/" + taus + "/C_" + str(ident) + ".dat", T)  # measured steady state density
        EV_comp[:, 0] = Ev.real
        EV_comp[:, 1] = Ev.imag
#
# 		for k in range(len(minpos)):
# 			target = []
# 			for m in range(MS[1]):
# 				target.extend(list(range(minpos[k] - rancut + m * MS[0] , minpos[k] + rancut+ m * MS[0])))  #potential dependent
# 			for t in range(len(minpos)):
# 				origin = []
# 				for m in range(MS[1]):
# 					origin.extend(list(range(minpos[t] - rancut + m*MS[0] , minpos[t] + rancut +m *MS[0])))
#
#
# 				tau = wline * 1. - 1.  # read from first line count...
# 				#mfpt_ar[t][k] = mfpt(L,target = target, origin = origin)
# 				############## DOES NOT WORK YET!!!!!!!!!!!! ##############
#
# 				# calculate the mean first passage time between(MFPT) origin and target
# 				# Origin and Target have to be specified, such that the represent the minimum of the used potential - general algorithm would be beneficial here. Could be a job for you.
#
# 		#MSM_ = MSM(L)
#
#
# 		length = 1000
# 		mfpt_ana = np.zeros((lvl*lvl,length))
#
#
# 		# calculatiion of first passage time distribution (FPTD) - Integrated FPTD should match MFPT
# 		# as before, start and end has to be specified just like before.
# 		for k in range(lvl):
# 			for t in range(lvl):
# 				if (k!=t):
# 					start = []
# 					for m in range(MS[1]):
# 						start.extend(list(range(minpos[k] - rancut + m * MS[0] , minpos[k] + rancut+ m * MS[0])))
# 					end = []
# 					for m in range(MS[1]):
# 						end.extend(list(range(minpos[t] - rancut + m * MS[0] , minpos[t] + rancut+ m * MS[0])))
#
# 					end = list(range(minpos[t]-rancut, minpos[t]+rancut+1))
# 					mfpt_ana[k*lvl+t,:] = first_passage(start,end,L,length)
#
# 		tau = wline - 1
# 		np.savetxt(path0+'single_particle/'+taus+'/MFPT/markovana_'+str(ident)+'.dat',np.transpose(mfpt_ana))
#
#
#
#
# 		MFPT_meas_reg = np.zeros((lvl,lvl))
# 		MFPT_err_reg = np.zeros((lvl,lvl))
# 		MFPT_steps_reg = np.zeros((lvl,lvl))
# 		MFPT_hist = np.zeros((lvl*lvl,1000))
#
# 		#print(tau, ident)
# 		#np.savetxt(path0+'single_particle/MFPT/markov_area_'+str(int(tau))+"_"+str(ident)+'.dat',(mfpt_ar))
#
#
# 		np.savetxt(path0+'single_particle/'+taus+'/T/T_'+str(ident)+'.dat',L)
# 		np.savetxt(path0+'single_particle/'+taus+'/EV/EV_com_'+str(ident)+'.dat',EV_comp)
# 		for k in range(Enum): #save complex eigenvectors up to order specified in the beginning of the script
# 			Evec_comp[:,0] = Evec[k,:].real
# 			Evec_comp[:,1] = Evec[k,:].imag
# 			for l in range(ms):
# 				if (Evec[0,l] > 0):
# 					Evecr_comp[l,0] = (Evec[k,l]/ Evec[0,l]).real
# 					Evecr_comp[l,1] = (Evec[k,l]/ Evec[0,l]).imag
##				else:
##					Evecr_comp[l,0] = 0.
##					Evecr_comp[l,1] = 0.
# 			np.savetxt(path0+'single_particle/'+taus+'/EV/Evec'+str(k+1)+'_com_'+str(ident)+'.dat',Evec_comp)
# 			np.savetxt(path0+'single_particle/'+taus+'/EV/Evecr'+str(k+1)+'_com_'+str(ident)+'.dat',Evecr_comp)
# 			# left and right eigenvectors only differ by state density - check in book
#
# 	for j in range(1,ms):
# 		if (Ev[j].real > 0 and not math.isclose(Ev[j].real,1.0)): #How to deal with complex part? < - this is an old comment from me - And I still dont know the answer.
# 			lt[j,i] = - tao / math.log(Ev[j].real)  # Only real part used here - so it actually wrong. So lagtime could be chosen badly, because the system is not markovian!
# 		else:
# 			lt[j,i] = 'NaN'
# 		EV_arr[j+1,i] = Ev[j].real #fill array with Ev's in order
#
# 	for k in range(Enum):
# 		Evec_arr[i,:,k] = Evec[k+1].real # Evec[0] cotains invariant distr.
# 		for l in range(ms):
# 			if (Evec[0,l].real > 0):
# 				Evecr_arr[i,l,k] = Evec[k+1,l].real / Evec[0,l].real #right Evec
# 			else:
# 				Evecr_arr[i,l,k] = 0.
#
#
# if db_flag: # do detailed balance test if wanted - detailed balance sould be satisfied in equilibrium
# 	db_mat = np.zeros((ms,ms))
# 	av_violation = 0.
# 	T_check = T_check.astype(float)
#
# 	T_check = T_check /  np.sum(T_check[:,:])	# normalise
#
# 	for i in range(ms):
# 		for j in range(ms):
# 			db_mat[i,j] = abs(T_check[i,j] -T_check[j,i])
# 			av_violation += db_mat[i,j]
# 	av_violation /= (ms*ms - ms)/2
#
# 	print("average violation of db: " ,av_violation)
#
# 	np.savetxt(path0+'single_particle/'+taus+'/detbal/err_'+str(ident)+'.dat',db_mat)
# 	np.savetxt(path0+'single_particle/'+taus+'/detbal/flux_'+str(ident)+'.dat',T_check)
#
#
# np.savetxt(path0+'single_particle/lagtime/'+str(ident)+'.dat',np.transpose(lt))
# np.savetxt(path0+'single_particle/EV/EV_'+str(ident)+'.dat',np.transpose(EV_arr))
#
#
#
# for k in range(Enum):
# 	np.savetxt(path0+'single_particle/EV/Evec'+str(k+1)+'_'+str(ident)+'.dat',np.transpose(Evec_arr[:,:,k]))
#
# for k in range(Enum):
# 	np.savetxt(path0+'single_particle/EV/Evecr'+str(k+1)+'_'+str(ident)+'.dat',np.transpose(Evecr_arr[:,:,k]))
#
#
#
#
#
#
#
#
#
#
#
