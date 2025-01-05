import argparse
import numpy as np
import math
import getpass


# classical ferrenberg swensen reweighting from ref potential

parser = argparse.ArgumentParser()
parser.add_argument('-o', help='ref')
parser.add_argument('-minT', help='minimum dT')
parser.add_argument('-maxT', help='maximum dT')
args = parser.parse_args();
ident = args.o
minT = float(args.minT)
maxT = float(args.maxT)
user = getpass.getuser()

minbeta = 1/maxT
maxbeta = 1/minT

if (user=="mbause"):
	path0 = "/u/mbause/data/"
	path1 = "/u/mbause/data/"
else:
	path0 = "/data/isilon/bause/"
	path1 = "/data/pckr194/bause/"


###########
potential = np.loadtxt(path0+"single_particle/potential_ms_"+str(ident)+".dat")
data = np.loadtxt(path0+"single_particle/100/EV/ps_"+str(ident)+".dat")
###############

potential = potential[:,1]
steps = 100
minpos = [1, 10, 18, 26] 
rc = 1
betastep = (maxbeta -minbeta) / steps

out  = np.zeros((steps,5))

for line in open(path0+"single_particle/param_"+str(ident)+".dat","r"):
	cont = line.split()
	if(cont[0] == 'T0'): # old version
		Temperature = float(cont[1])
	if(cont[0] == 'T'): # new version
		Temperature = float(cont[1])


betaref = 1./Temperature

for i in range(steps):
	beta = minbeta+ i * betastep	

	rewdata = data * np.exp(-(beta-betaref) * potential)	

	rewdata = rewdata / (sum(rewdata))
	out[i,0] = 1./beta
	for j in range(4):
		out[i,j+1] = sum(rewdata[minpos[j]-rc:minpos[j]+rc+1])

np.savetxt(path0+"/single_particle/rew_"+ident+".dat",out)

#####################################
P = np.zeros((10,5))
ident = 8050
for i in range(10):
	ide = (int(ident) + i*100)
	data = np.loadtxt(path0+"single_particle/100/EV/ps_"+str(ide)+".dat")

	for line in open(path0+"single_particle/param_"+str(ide)+".dat","r"):
		cont = line.split()
		if(cont[0] == 'T0'): # old version
			Temperature = float(cont[1])
		if(cont[0] == 'T'): # new version
			Temperature = float(cont[1])

	P[i,0] = Temperature
	for j in range(4):
		P[i,j+1] = sum(data[minpos[j]-rc:minpos[j]+rc+1])


print(P[9,:])	
print(out[0,:])	

np.savetxt(path0+"/single_particle/rew_meas_"+str(ident)+".dat",P)













