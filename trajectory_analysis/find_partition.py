import numpy as np
import math
import argparse
import os
from define_flow import define_flow_direction
from define_flow import read_cut
from define_flow import read_minpos
from define_flow import pos_pot

def import_timeseries(path):
	#use coursed timeseries for now
	a = []
	v = []
	for line in open(path, "r"):
		stuff = line.split()
		a.append(float(stuff[0]))
		v.append(float(stuff[1]))

	a = np.asarray(a)
	v = np.asarray(v)
	return a,v
	

	
def split_series(series, T = 10):
	# T is number of steps
	N_tr = math.floor(len(series) / T)
	Gamma = series.reshape((N_tr,T))	
	return Gamma

def get_pos(x,ms):
	pos = math.floor(x*ms)
	if(pos == ms):
		 pos -= 1
	return pos

def J(tra,v,ms,T):
	flux = 0
	pos1 = get_pos(tra[0],ms)	
	for i in range(1,T+1):
		pos2 = get_pos(tra[i],ms)	
		if (pos2 != pos1):
			if (v[i] > 0):
				if(pos2 > pos1):
					flux += (pos2-pos1)	
				else:
					flux += ms - pos1 + pos2 #b.c.
			else:
				if(pos1 > pos2):
					flux += (pos2-pos1) # becomes negative
				else:
					flux += pos2 - pos1 -ms #b.c.
			pos1 = pos2

	return flux

def Jloc(pos1,pos2,v1,ms):
	flux = 0
	if (pos2 != pos1):
		if (v1 > 0):
			if(pos2 > pos1):
				flux =1
			else:
				flux = -1
		else:
			if(pos1 > pos2):
				flux = -1
			else:
				flux = 1

	return flux


def get_E(tra,pot,ms,T):
	En = 0
	for i in range(1,T+1):
		posn = get_pos(tra[i],ms)
		En += pot[pos] # table with E values is needed here
	
	En /= T # scale to average Energy
	return En


def get_J(Gamma,vGamma, ms):
	N = len(Gamma[:,0])
	J_vec = np.zeros(N)
	for i in range(N):
		J_vec[i] = J(Gamma[i,:],vGamma[i,:],ms)			
		
	return J_vec

def create_hist(J):
	N = len(J)
	jmax = int(round(max(J)))
	jmin = - int(round(min(J)))
	hist = np.zeros(jmax + jmin + 1)	
	histx = np.arange(-jmin,jmax+1)	
	for i in range(len(J)):
		hist[J[i] + jmin] = hist[J[i] + jmin] +1
	
	return histx,hist

#def iterate_Z(hist,Nk):
#	num = len(hist[:,0])
#	Z = np.ones(num) # init
#	Zold = np.ones(num) # init
#	cut = 1
#	eps = 0.1
#	while (cut > eps): #iteration loop
#		for i in range(num): #loop over all hists
#			for n in range(Nk[i]): #loop over all J's
#				for k in range(num):
				
def get_stuff(Gamma,v,pot, ms,T):
	flux = J(Gamma,v,ms,T)	
	Energy = get_E(Gamma,pot,ms,T)		
	#for i in range(len(beta)):
	#	for j in range(len(lam)):
			#out[i,j] = math.exp(-lam[j] * flux - beta[i] * Energy)
			
	return flux, Energy

def define_cross(minpos,dist,ms): 
	states = len(minpos)
	cross = np.zeros((states,2))
	if (cut[0] <  minpos[0]):
		for k in range(2,states):
			cross[k,0] = min(minpos[k-2]+dist+1,minpos[k]-dist)
			cross[k,1] = max(minpos[k-2]+dist,minpos[k]-dist-1)
		
		cross[0,1] = min(minpos[0]-dist-1,minpos[1]+dist)
		cross[0,0] = max(minpos[0]-dist,minpos[1]+dist+1)
		cross[1,1] = min(minpos[1]-dist-1,minpos[states-1]+dist)
		cross[1,0] = max(minpos[1]-dist,minpos[states-1]+dist+1)
		
		outpos = np.where(cross >= ms)
		cross[outpos] -= ms

		outpos = np.where(cross < 0)
		cross[outpos] += ms

	else:	
		print("TODO : case cut[0] > minpos[0] ")			
	
	print(cross)

	return cross	

#better trajectory data needed.
#T = 10



parser = argparse.ArgumentParser()
parser.add_argument('-o','--output', help='ident of output')
parser.add_argument('-N','--number', help='number of input data')
args = parser.parse_args();
ident = args.output
Nmeas = args.number

ms,minpos = read_minpos(ident)
ms,cut = read_cut(ident)


print("min",minpos)
print("cut",cut)

blocks = 10

MFPT = np.zeros((blocks,len(cut),len(cut)))
MFPT_steps = np.zeros((blocks,len(cut),len(cut)))


maxstep = 20000
MFPT_hist = np.zeros((9,maxstep))

path = "/home/theorie/bause/code/single_particle/data/potential_" +str(ident) +".dat"
x,pot = np.loadtxt(path, usecols=(0,1),unpack=True, delimiter="\t")

pot_v = np.zeros(ms)
pos_old = get_pos(x[0],ms)
U = 0
count =0
for i in range(1, len(x) ):
	pos = get_pos(x[i],ms)
	if(pos  == pos_old):
		U = (U+pot[i]) / (count +1)
		count +=1
		pot_v[pos] = U
	else:
		count =1
		U = pot[i]
	pos_old = pos


path = "/data/pckr146/bause/single_particle/trajectory_" +str(ident) +".dat"
#data,v = import_timeseries(path)

series = open(path, "r")

datalen = int(os.popen("wc -l "+path).readline().split()[0])

print(datalen)
datalen = int(Nmeas)

print(datalen)

datalen_block = int(datalen /blocks)

histN= ms
maxT = 500
maxit = 10
Tstep = int(maxT / maxit)



num_lambda = 12
num_beta = 10
delta_lambda = 0.1
delta_beta = 0.1
Q = np.zeros((num_beta, num_lambda))

E_hist = np.zeros(ms)
flux_hist= np.zeros(3)

lag_lambda = list(delta_lambda * k for k in range(num_lambda))
lag_beta = list(delta_beta * k for k in range(num_beta))


Q = np.ones((maxit, num_beta, num_lambda)) #trajectory of length 1 is useles
Gamma = np.zeros(maxT+1) # one more for old one
v_part = np.zeros(maxT+1)
pos_old = ms /2


ll =2
long_hist=np.zeros((maxit,2*ll*maxT+1))

for i in range(1000): #start in the middle
    series.readline()
	
stuff = series.readline()
split = stuff.split()
Gamma[0] = float(split[0])
v_part[0] = float(split[0])
pos = get_pos(Gamma[0],ms)
E_hist[pos] +=1

pos_old = pos


	
for m in range(1, maxT+1): #initial trajectory
	stuff = series.readline()
	split = stuff.split()
	Gamma[m] = float(split[0])
	v_part[m] = float(split[1])
	pos = get_pos(Gamma[m],ms)
	E_hist[pos] +=1
	fl = Jloc(pos,pos_old,v_part[m],ms) +1
	flux_hist[fl] += 1
	pos_old = pos
	


#fl = J(Gamma,v_part,ms,maxT)
lvl = len(cut)

#long_hist[fl+maxT] +=1
# change such that several histograms for different tra-length are recorded
# does it change anything?

rancut = 0 #pretty random

cross = define_cross(minpos,rancut,ms)

reg1 = pos_pot(pos_old,lvl,cut)
if (reg1 == lvl):
	reg =0
else:
	reg = reg1

steps = np.zeros((lvl,lvl))



reg = pos_pot(pos,lvl,cut)  # MFPT stuff

if(reg == lvl):
	reg = 0

regmax = lvl -1

l = 10
for b in range(0,blocks):
	for j in range(0,datalen_block):
		if (abs(pos_old - pos) > ms/2):
			if(pos < pos_old):
				ran1 = range(-2,pos+1)
				ran2 = range(pos_old,ms+2)
			else:
				ran1 = range(-2,pos_old+1)
				ran2 = range(pos,ms+2)
		else:
				if (pos > pos_old):
					ran1 = range(pos_old,pos+1)
					ran2 = ran1
				else:
					ran1 = range(pos,pos_old+1)
					ran2 = ran1
			

		if (cross[reg,1] in ran1 or cross[reg,1] in ran2 or cross[reg,0] in ran1 or cross[reg,0] in ran2 ):
		# if new region was reached and minimum position was crossed	

			reg_old = reg
	
			if(cross[reg,0] in ran1 or cross[reg,0] in ran2):
				reg =  reg - 1 if reg > 0 else regmax
			else:
				reg = reg +1 if reg < regmax else 0

			#print(b*datalen_block+j,reg_old, reg, pos_old, pos, steps, Gamma[maxT])

			MFPT[b,reg_old,reg] = (MFPT[b,reg_old,reg] * MFPT_steps[b,reg_old,reg] + steps[reg_old,reg]) / (MFPT_steps[b,reg_old,reg] +1.)
			MFPT_steps[b,reg_old,reg] +=1 
			if (steps[reg_old,reg] > maxstep):
				print("short hist")
			else:
				MFPT_hist[reg_old*3 + reg][np.int(steps[reg_old,reg])] += 1
			steps[:,:] +=1
			steps[reg_old,reg] = 1
		else:
			steps[:,:] += 1

		for i in range(2,maxit+2):
			tr_length = int((i-1) * Tstep)
			#print(j,tr_length)
			if (j % tr_length == 0):
				fl, energy = get_stuff(Gamma,v_part,pot_v,ms,tr_length)	
				#fl = J(Gamma,v_part,ms,tr_length)
				if (abs(fl) > maxT*ll):
					#fl = tr_length
					print("problem" ,Gamma)
					print(v_part,fl,i)
				else:
					long_hist[i-2,fl+ll * maxT] +=1 

	
				for k in range(num_lambda):
					for l in range(num_beta):				
						Q[i-2,l,k] = np.logaddexp(Q[i-2,l,k], -lag_beta[l] * energy - lag_lambda[k] * fl  / tr_length )

	#if (j % 10000 ==0):
	#	print(j)

		Gamma=np.delete(Gamma,0,0)
		v_part=np.delete(v_part,0,0)
		stuff = series.readline()
		split = stuff.split()
		Gamma=np.append(Gamma,float(split[0]))
		v_part=np.append(v_part,float(split[1]))

		pos_old = pos
		pos = get_pos(Gamma[maxT],ms)
		E_hist[pos] +=1
		fl = Jloc(pos_old,pos,v_part[maxT],ms) +1
		flux_hist[fl] += 1	



for i in range(maxit-1):
	np.savetxt("/data/isilon/bause/single_particle/partition_fct/"+str(ident)+"_"+str(i)+".dat",Q[i])

print(np.shape(MFPT_hist))
for i in range(0,9):
	denom = sum(MFPT_hist[i,:])
	if (denom > 0 ):
		MFPT_hist[i,:] = MFPT_hist[i] / denom	

print(np.shape(MFPT_hist))

MFPT_err =np.zeros((len(cut),len(cut)))
MFPT_av = np.zeros((len(cut),len(cut)))
MFPT_s = np.zeros((len(cut),len(cut)))
for b in range(blocks):
	MFPT_av[:,:] += MFPT[b,:,:]	
	MFPT_s += MFPT_steps[b,:,:]

MFPT_av = MFPT_av / blocks



for b in range(blocks):
	fac = MFPT[b,:,:] - MFPT_av[:,:] 
	MFPT_err += fac* fac

MFPT_err = np.sqrt( MFPT_err/ (blocks *(blocks-1)))



	
np.savetxt("/data/isilon/bause/single_particle/MFPT/hist_"+str(ident)+".dat",np.transpose(MFPT_hist))

np.savetxt("/data/isilon/bause/single_particle/MFPT/meas_"+str(ident)+".dat",MFPT_av)
np.savetxt("/data/isilon/bause/single_particle/MFPT/measerr_"+str(ident)+".dat",MFPT_err)
np.savetxt("/data/isilon/bause/single_particle/MFPT/measnum_"+str(ident)+".dat",MFPT_s)
np.savetxt("/data/isilon/bause/single_particle/partition_fct/Ehist_"+str(ident)+".dat",E_hist)
np.savetxt("/data/isilon/bause/single_particle/partition_fct/fluxhist_"+str(ident)+".dat",flux_hist)

# make sure that flux is first line in long_hist
Js = np.arange(2*ll*maxT+1)
Js = Js - ll*maxT
out = np.concatenate((np.asarray([Js]), np.asarray(long_hist)),axis = 0)

firstline = [str(i*Tstep)  for i in range(0,maxit+1)]
firstline[0] = "#" 

np.savetxt("/data/isilon/bause/single_particle/partition_fct/fluxlong_"+str(ident)+".dat",np.transpose(out))

with open("/data/isilon/bause/single_particle/partition_fct/fluxlong_"+str(ident)+".dat", "r+") as f:
	f.seek(0)
	for item in firstline:	
		f.write("%s\t" % item)
	f.write("\n")



#for i in range(num_in):
#	num = 4000 + i * 100 #more elegant
#	path = "/data/isilon/bause/single_particle/trajectory_" +str(num) +".dat"
#	data,v = import_timeseries(path)
#	Gamma = split_series(data,T)
#	vGamma = split_series(v,T)

#	J = get_J(Gamma,vGamma,ms)
#	Jx,p = create_hist(J)
#	Nk = sum(p)
#	p = p/Nk
#	Jmin = min(Jx,Jmin)
	# order hist such that histograms are at same pos
#	hist[i] = p

#np.savetxt("test2.dat", hist)













