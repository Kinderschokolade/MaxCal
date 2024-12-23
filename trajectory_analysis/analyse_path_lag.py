import numpy as np
import os
import math
import argparse

def H(x, kc,shift):
	return(1./(1.+np.exp(-2*kc*(x-shift))));



def dH(x, kc,shift):
	return(-2.*np.exp(-2.*kc*(-shift+x))*kc)/((1.+np.exp(-2*kc*(-shift+x)))*(1.+np.exp(-2*kc*(-shift+x))));


def F_tanh(x, kc, U, xc):
	lvl = np.size(U)
	f = (dH(x,kc,xc[0]) - dH(x,kc,xc[1])+ dH(x,kc,xc[lvl]))* U[0] +(dH(x,kc,xc[lvl-1]) - dH(x,kc,xc[lvl]) -dH(x,kc,xc[0]))* U[lvl-1]

	for i in range(1,lvl-1):
		f += (dH(x,kc,xc[i]) - dH(x,kc,xc[i+1]))* U[i]
	return f

def get_U(x, kc, U, xc):
	E = 0.
	lvl = len(U);

	E = (H(x,kc,xc[0]) - H(x,kc,xc[1])+ H(x,kc,xc[lvl]))* U[0] +(H(x,kc,xc[lvl-1]) - H(x,kc,xc[lvl]) - H(x,kc,xc[0]))* U[lvl-1] ;
	
	for j in range(1,lvl-1):
		E += (H(x,kc,xc[j]) - H(x,kc,xc[j+1]))* U[j] ;
		
	return E

def get_pos(x,ms):
	pos = math.floor(x*ms)
	if(pos == ms):
		 pos -= 1 # rare case
	return pos

def get_core_pos(x,ms,shift):
	pos = math.floor((x+shift)*ms*2)
	if(pos%2 == 0):
		pos = pos//2
	else:
		pos = ms
	return pos


parser = argparse.ArgumentParser()
parser.add_argument('-o', help='ident')
parser.add_argument('-ms', help='end', default=30)
args = parser.parse_args();
ident = args.o
ms = int(args.ms)




path = "/data/pckr194/bause/single_particle/trajectory_" +str(ident) +".dat"
series = open(path, "r")
tra_len = int(os.popen("wc -l "+path).readline().split()[0])

U = np.zeros(6)
xc = np.zeros(7)

for line in open("/data/isilon/bause/single_particle/param_"+str(ident)+".dat","r"):
	cont = line.split()
	if(cont[0] == 'dT'):
		dT = float(cont[1])
	if(cont[0] == 'gamma'):
		gamma = float(cont[1])
	if(cont[0] == 'T0'):
		Temperature = float(cont[1])
	if(cont[0] == 'extf'):
		force = float(cont[1])
	if(cont[0] == 'U0'):
		U[0] = float(cont[1])
	if(cont[0] == 'U1'):
		U[1] = float(cont[1])
	if(cont[0] == 'U2'):
		U[2] = float(cont[1])
	if(cont[0] == 'U3'):
		U[3] = float(cont[1])
	if(cont[0] == 'U4'):
		U[4] = float(cont[1])
	if(cont[0] == 'U5'):
		U[5] = float(cont[1])
	if(cont[0] == 'xc0'):
		xc[0] = float(cont[1])
	if(cont[0] == 'xc1'):
		xc[1] = float(cont[1])
	if(cont[0] == 'xc2'):
		xc[2] = float(cont[1])
	if(cont[0] == 'xc3'):
		xc[3] = float(cont[1])
	if(cont[0] == 'xc4'):
		xc[4] = float(cont[1])
	if(cont[0] == 'xc5'):
		xc[5] = float(cont[1])
	if(cont[0] == 'xc6'):
		xc[6] = float(cont[1])
	if(cont[0] == 'k'):
		kc = float(cont[1])

shift = -1./(ms*4) # 50 % core position 
stuff =series.readline()
split = stuff.split()
x = float(split[0])
pos = get_core_pos(x,ms,shift)
#pos = get_pos(x,ms)
pos_old = pos
Sprod = 0.
Sprod2 = 0.
SprodA = 0.



sqdT = np.sqrt(dT)
sigma = np.sqrt(2*Temperature/gamma)



traj = [] 
ran_l =[]
F_l = []
pos_l = []
V_l = []

checkpos = np.zeros((ms,ms))
steps = np.zeros((ms,ms))

lag = 80
count =0
while (count < lag):
	stuff = series.readline()
	if stuff.strip():
		split = stuff.split()
		xold = x
		x = float(split[0])
		ran = float(split[2])
		F= F_tanh(x,kc,U,xc) + force
		V= get_U(x,kc,U,xc)
		pos = get_core_pos(x,ms,shift)
		#pos = get_pos(x,ms)
		steps = steps +1
		if(pos != pos_old and pos != ms):
			pos_old = pos
			checkpos[pos,:] = 1
			steps[pos,:] = 1  
			for l in range(ms):
				if(checkpos[l,pos] ==1):
					steps[l,pos] =0
					checkpos[l,pos] = 0

		traj.append(x)
		ran_l.append(ran)
		F_l.append(F)		
		pos_l.append(pos)
		V_l.append(V)
		count = count +1

	
path1 = "/data/pckr194/bause/traj/tr1_"
path2 = "/data/pckr194/bause/traj/tr2_"

Sprod_th = np.zeros((ms,ms))
for i in range(ms):
	for j in range(ms):
		Sprod_th[i,j] = ( -force / ms * (i-j) + get_U(i/ms,kc,U,xc) - get_U(j/ms,kc,U,xc)) / Temperature
		print(i,j,Sprod_th[i,j])
Sprod_data= []
Action_data= []
for i in range(0,ms):
    Sprod_data.append([])
    Action_data.append([])
    for j in range(0,ms):
        Sprod_data[i].append([])
        Action_data[i].append([])

#tlen = 0

steps.astype(int)

fac1 = sqdT/(gamma*sigma)
fac2 = dT/(2*gamma*gamma*sigma*sigma)

count =0
ccount =0


tracount = int(tra_len/10);
Action_cur = 1000

for i in range(tra_len):
	stuff = series.readline()
	if stuff.strip():
		if (i % tracount ==0):
			print(i / tracount)
		split = stuff.split()
		x = float(split[0])
		ran = float(split[2])
		F= F_tanh(x,kc,U,xc) + force
		V= get_U(x,kc,U,xc)
		pos = get_core_pos(x,ms,shift)
		#pos = get_core_pos(x,ms)
		steps = steps+1 #for all possible transitions		
		traj.append(x)
		ran_l.append(ran)
		F_l.append(F)
		pos_l.append(pos)
		V_l.append(V)
		Sprod = 0 
		Action = 0

		l = pos_l[0]
		for p in range(1,lag): # length-1 steps
			#F= F_l[j]+ F_l[lag - j+1]
			po = (traj[p-1]+traj[p])/2.
			Action += (traj[p] - traj[p-1])*(traj[p]-traj[p-1])/(dT*2.) + V_l[p]*dT  # nothing about f_ext yet....
			Sprod += (traj[p] - traj[p-1])*(F_tanh(po,kc,U,xc))/(Temperature)  # stratonovic



		if (l ==17 and pos ==19): # and Action < Action_cur):
			count = count+1
			#Action_cur = Action
			print(count, l,pos, Sprod,Action,Sprod_th[17][19])
			np.savetxt("/data/pckr194/bause/traj/tr1_"+str(count)+".dat",traj)

		if (l ==19 and pos ==17): # and Action < Action_cur):
			ccount = ccount+1
			#Action_cur = Action
			print(ccount, l,pos, Sprod,Action,Sprod_th[19][17])
			np.savetxt("/data/pckr194/bause/traj/tr2_"+str(ccount)+".dat",traj)



		#Sprod_data[pos_l[0]][pos].append(Sprod)
		#Action_data[pos_l[0]][pos].append(Action)
			
		del traj[0]
		del ran_l[0]
		del F_l[0]
		del pos_l[0]
		del V_l[0]


#print("min: ", pos_min , Sprod_min, Sprod_th)
#print("max: ", pos_max , Sprod_max, Sprod_th)
#
#
#bins = 200
#hist= np.zeros((bins,2))
#for i in range(ms):
#	for j in range(ms):
#		if Sprod_data[i][j]:
#			ra = 4* abs(Sprod_th[i,j])
#			hist[:,1],edgesf = np.histogram(np.asarray(Sprod_data[i][j]),bins = bins)
#			for k in range(bins):
#				hist[k,0] = (edgesf[k] + edgesf[k+1]) / 2.
#			np.savetxt("/data/isilon/bause/single_particle/Sprod/hist_"+str(i)+"_"+str(j)+"_"+str(ident)+".dat", hist)



