import numpy as np
import argparse
import math
import os
from analyse_trajectory import mfpt_trajectory_cts_cross
from analyse_trajectory import get_pos
from analyse_trajectory import get_reg

from define_flow import read_cut
from define_flow import read_minpos

def flow_direction(pl, mi, ms):
	r = np.zeros((ms,ms))
	for i in range(ms):
		for j in range(ms):
			if (pl[i,j] > mi[i,j]):
				r[i,j] = +1
			else:
				r[i,j] = -1

	return r


parser = argparse.ArgumentParser()
parser.add_argument('-o', help='ident')
parser.add_argument('-N', help='trajectory length')
parser.add_argument('-l', help='lag')
parser.add_argument('-rc', help='rancut')
args = parser.parse_args();
ident = args.o
Nmeas = args.N
lag = int(args.l)
rancut = int(args.rc)



ms = 0
dT = 0
for line in open("/data/isilon/bause/single_particle/param_"+str(ident)+".dat","r"):
	cont = line.split()
	if(cont[0] == 'microstates'):
		ms = int(cont[1])
	if(cont[0] == 'T0'):
		T = float(cont[1])	
	if(cont[0] == 'dT'):
		dT = float(cont[1])
	if(cont[0] == 'extf'):
		f = float(cont[1])
	

path = "/data/isilon/bause/"
ms, minpos = read_minpos(path,ident)
ms, cut = read_cut(path,ident)

lvl = len(minpos)
blocks = 10
#print(minpos)

steps = np.zeros((lag,ms,ms))


maxstep =10000
FPT_hist = np.zeros((lag,lvl*lvl,maxstep))  # maximal steps is 1000....
reg_steps = np.zeros((lvl,lvl))


path = "/data/pckr194/bause/single_particle/trajectory_" +str(ident) +".dat"

series = open(path, "r")

datalen = int(os.popen("wc -l "+path).readline().split()[0])

therm_len = 1000

datalen = int(Nmeas)
datalen = int((datalen - 2*therm_len ) / lag)
#print(datalen)

datalen_block = int(datalen /blocks)



for i in range(therm_len): #ignore first therm steps
    series.readline()


checkpos = np.zeros((ms,ms))
checkreg = np.zeros((lvl,lvl))

stuff = series.readline()
split = stuff.split()
x = float(split[0])
v = float(split[1])
pos = get_pos(x,ms)
reg = get_reg(pos, minpos, rancut)
pos_old = pos
reg_old = reg
checkpos[pos,:] = 1
checkreg[reg,:] = 1
xold = x

Sp = np.zeros(blocks)

for s in range(0,blocks):
	for j in range(0,datalen_block):
#		#print(x, pos_old, pos, v)
		if (pos != pos_old):	
			reg = get_reg(pos, minpos, rancut)
			if reg is not None:
				checkreg[reg,:] = 1

			for rego_it in range(lvl):
				if reg is not None:
					if (checkreg[rego_it,reg] ==1 and rego_it != reg):
						if (reg_steps[rego_it,reg]<maxstep):
							FPT_hist[s,rego_it*lvl+reg,int(reg_steps[rego_it,reg])] +=1
						else:
							print("dgapkjbs")
						reg_steps[rego_it,reg] = 0 # reset counter
						checkreg[rego_it,reg] = 0 # trajectory reg_old -> reg ist now inactive
						checkreg[reg,:] = 1 # start new trajectory reg -> all other reg

		for i in range(ms):
			for j in range(ms):
				if (checkpos[i,j] == 1):
					steps[i,j] += 1
		for i in range(lvl):
			for j in range(lvl):
				if (checkreg[i,j] == 1):
					reg_steps[i,j] += 1

		Sp[s]+= f * v /T

		#update to next step according to lag
		v =0
		for i in range(lag):
			stuff = series.readline()
			split = stuff.split()
			v += float(split[1])

		xold = x
		x = float(split[0])

		pos_old = pos
		pos = get_pos(x,ms)
		reg = get_reg(pos,minpos,rancut)

#	print(MFPT[b,0,1], MFPT_steps[b,0,1])

#MFPT_av = np.zeros((ms,ms))
#MFPT_err = np.zeros((ms,ms))
#MFPT_stepsav = np.zeros((ms,ms))

#for b in range(blocks):
#	MFPT_av += MFPT[b,:,:]
#	MFPT_stepsav += MFPT_steps[b,:,:]

#MFPT_av /= blocks

#for b in range(blocks):
#	temp = MFPT_av - MFPT[b,:,:]
#	MFPT_err += temp *temp
#MFPT_err = np.sqrt(1/(blocks*(blocks-1))* MFPT_err  )

Sp = Sp / datalen_block
Sp_av = sum(Sp) / blocks

Sp_sig = sum((Sp_av - Sp[i])*(Sp_av - Sp[i]) / blocks for i in range(blocks))

print(Sp_av, np.sqrt(Sp_sig))

MFPT_cts_av = np.zeros((lvl,lvl))
MFPT_cts_err = np.zeros((lvl,lvl))

FPT_av = np.zeros((lvl*lvl,maxstep))
FPT_err = np.zeros((lvl*lvl,maxstep))


FPT_av = sum(FPT_hist[i,:,:] for i in range(lag))
FPT_err = sum(FPT_hist[i,:,:] *FPT_hist[i,:,:] for i in range(lag))
#print(FPT_err)
FPT_err = np.sqrt( FPT_err/lag - FPT_av *FPT_av/ (lag*lag))

#print(FPT_err)



for i in range(lvl):
	for j in range(lvl):
		if (i != j):
			denom = sum(FPT_av[i*lvl+j,:])
			FPT_av[i*lvl+j,:] = FPT_av[i*lvl+j,:] / denom
			FPT_err[i*lvl+j,:] = FPT_err[i*lvl+j,:] / denom

for i in range(lvl):
	for j in range(lvl):
		if (i != j):
			MFPT_cts_av[i,j] = sum( FPT_av[i*lvl+j,k] * k for k in range(maxstep))	
			MFPT_cts_err[i,j] = sum( FPT_av[i*lvl+j,k] * k *k for k in range(maxstep))	

MFPT_cts_err = np.sqrt(MFPT_cts_err - MFPT_cts_av * MFPT_cts_av)

#print(MFPT_cts_av)

tau_s = str(int(lag))

#np.savetxt("/data/isilon/bause/single_particle/"+tau_s+"/MFPT/ctsmeas_"+str(ident)+".dat",MFPT_av)
#np.savetxt("/data/isilon/bause/single_particle/"+tau_s+"/MFPT/ctsmeaserr_"+str(ident)+".dat",MFPT_err)

path = "/data/isilon/bause/single_particle/"
temppath = "/home/theorie/bause/data/"

np.savetxt(temppath+tau_s+"/MFPT/ctsmeas_"+str(ident)+".dat",MFPT_cts_av)
np.savetxt(temppath+tau_s+"/MFPT/ctsmeaserr_"+str(ident)+".dat",MFPT_cts_err)

np.savetxt(temppath+tau_s+"/MFPT/ctshist_"+str(ident)+".dat",np.transpose(np.concatenate((FPT_av,FPT_err))))




