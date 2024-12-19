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

def create_hist(data,bins,minent):

	hist = np.zeros((bins,2))
	hist[:,1],edgesf = np.histogram(np.asarray(data),bins = bins, density = False)
	hist[:,1] = hist[:,1] / sum( hist[:,1] )
	pos = np.squeeze(np.where(hist[:,1]>minent))

	if(pos.size > 1):
		ran = (edgesf[pos[0]], edgesf[pos[-1]+1])
		hist[:,1],edgesf = np.histogram(np.asarray(data),bins = bins,range = ran, density=False )
		hist[:,1] = hist[:,1] / sum( hist[:,1] )

	for k in range(bins):
		hist[k,0] = (edgesf[k] + edgesf[k+1]) / 2.

	return hist
	


parser = argparse.ArgumentParser()
parser.add_argument('-o', help='ident')
parser.add_argument('-ms', help='end', default=30)
parser.add_argument('-l','--lag', nargs='+', help='lagtimes', required=True)
args = parser.parse_args();
ident = args.o
ms = int(args.ms)
lagv = args.lag


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
pos_old = pos

print(Temperature)


sqdT = np.sqrt(dT)
sigma = np.sqrt(2*Temperature/gamma)



traj = [] 
ran_l =[]
F_l = []
Fm_l = []
pos_l = []
V_l = []


lag = [int(lagv[i]) for i in range((len(lagv)))] 
lagc = len(lag)
lagmax = max(lag)
count =0
lag.extend([0])
print("lagcheck", lag)
while (count < lagmax):
	stuff = series.readline()
	if stuff.strip():
		split = stuff.split()
		xold = x
		x = float(split[0])
		ran = float(split[2])
		F= F_tanh(x,kc,U,xc) + force
		Fm= F_tanh((xold+x)/2.,kc,U,xc) + force
		V= get_U(x,kc,U,xc)
		pos = get_pos(x,ms)
		pos_l.append(pos)
		traj.append(x)
		ran_l.append(ran)
		F_l.append(F)		
		Fm_l.append(F)		
		V_l.append(F)		
		count = count +1
	
path1 = "/data/pckr194/bause/traj/tr1_"
path2 = "/data/pckr194/bause/traj/tr2_"

Sprod_th = np.zeros((ms,ms))
delta = 1/(2*ms)
for i in range(ms):
	for j in range(ms):
		Sprod_th[i,j] = ( -force / ms * (i-j) + get_U(i/ms+delta,kc,U,xc) - get_U(j/ms+delta,kc,U,xc)) / Temperature
		print(i,j,Sprod_th[i,j])

Sprod_data= []
Action_data= []
SprodI_data= []
SprodS_data= []
SprodS2_data= []
SprodG_data= []
SprodG2_data= []
SprodG1_data= []
for i in range(0,lagc):
	Sprod_data.append([])
	Action_data.append([])
	SprodI_data.append([])
	SprodS_data.append([])
	SprodS2_data.append([])
	SprodG2_data.append([])
	SprodG1_data.append([])
	SprodG_data.append([])
	for j in range(0,ms):
		Sprod_data[i].append([])
		Action_data[i].append([])
		SprodI_data[i].append([])
		SprodS_data[i].append([])
		SprodS2_data[i].append([])
		SprodG_data[i].append([])
		SprodG1_data[i].append([])
		SprodG2_data[i].append([])
		for k in range(0,ms):
			Sprod_data[i][j].append([])
			Action_data[i][j].append([])
			SprodI_data[i][j].append([])
			SprodS_data[i][j].append([])
			SprodS2_data[i][j].append([])
			SprodG_data[i][j].append([])
			SprodG1_data[i][j].append([])
			SprodG2_data[i][j].append([])
	
fac2 = dT/(gamma*Temperature)
fac3 = np.sqrt(8.*dT/(gamma*Temperature))
fac4 = 2.
fac5 = sqdT/(gamma*sigma)
fac6 = dT/ (2.*gamma*gamma*sigma*sigma)

fac1 = np.sqrt(2.*Temperature*gamma/dT) 

count =0


tracount = int(tra_len/10);

for i in range(tra_len):
	stuff = series.readline()
	if stuff.strip():
		if (i % tracount ==0):
			print(i/ tracount)
		split = stuff.split()
		x = float(split[0])
		ran = float(split[2])
		F= F_tanh(x,kc,U,xc) + force
		V= get_U(x,kc,U,xc)
		pos = get_pos(x,ms)
		traj.append(x)
		ran_l.append(ran)
		F_l.append(F)
		F= F_tanh((x+traj[lagmax-1])/2.,kc,U,xc) + force
		Fm_l.append(F)  
		V_l.append(F)
		pos_l.append(pos)
		pos_old = pos_l[0]
		Sprod = 0 
		Action = 0
		SprodI = 0
		SprodS = 0
		SprodS2 = 0
	
		for k in range(lagc): # Girsanaov requires to repeat each tra from the beginning
			SprodG = 0
			SprodG2 = 0
			for p in range(lag[k]): # length-1 steps	
				F= +F_l[lag[k]-p] - F_l[p] #- 2. * force  Something does not work with the Girsanov stuff
				SprodG= -F *fac5 *ran_l[p]
				SprodG2+= F*F *fac6

			for p in range(lag[k-1]+1, lag[k]+1):
				#Sprod1 += fac2* F_l[p] * F_l[p]
				#Sprod2 += fac3* F_l[p] * ran_l[p]
				#Sprod3 += fac4* ran_l[p] * ran_l[p]
				Sprod +=  fac2* F_l[p] * F_l[p] # + fac3* F_l[p] * ran_l[p] + fac4* ran_l[p] * ran_l[p]
				Action += (traj[p] - traj[p-1])*(traj[p]-traj[p-1])/(dT*2.) + V_l[p]*dT  # nothing about f_ext yet....
#				SprodA2 += (traj[p] - traj[p-1])*(F_l[p-1]+F_l[p]+fac1*(ran_l[p]+ran_l[p-1])  )/(2.*Temperature)
#				SprodA += (traj[p] - traj[p-1])*(F_l[p-1]+ fac1*ran_l[p]  )/(Temperature)
				SprodS += (traj[p] - traj[p-1])*(Fm_l[p-1])/(Temperature)  # stratonovic
				SprodS2 += (traj[p] - traj[p-1])*(F_l[p-1]+F_l[p])/(2.*Temperature)  # stratonovic
				SprodI += (traj[p] - traj[p-1])*(F_l[p-1] )/(Temperature)

			Sprod_data[k][pos_old][pos_l[lag[k]]].append(Sprod)
			Action_data[k][pos_old][pos_l[lag[k]]].append(Action)   
			SprodS_data[k][pos_old][pos_l[lag[k]]].append(SprodS)   
			SprodS2_data[k][pos_old][pos_l[lag[k]]].append(SprodS2)   
			SprodI_data[k][pos_old][pos_l[lag[k]]].append(SprodI)   
			SprodG_data[k][pos_old][pos_l[lag[k]]].append(SprodG+SprodG2)   
			SprodG1_data[k][pos_old][pos_l[lag[k]]].append(SprodG)   
			SprodG2_data[k][pos_old][pos_l[lag[k]]].append(SprodG2)   
	
			
		del traj[0]
		del ran_l[0]
		del F_l[0]
		del pos_l[0]
		del V_l[0] 
		del Fm_l[0] 



#print("min: ", pos_min , Sprod_min, Sprod_th)
#print("max: ", pos_max , Sprod_max, Sprod_th)
#
#

bins = 100
minent = 1./(bins*3)


for k in range(lagc):
	for i in range(ms):
		for j in range(ms):
			if Sprod_data[k][i][j]:
				hist  = create_hist(Sprod_data[k][i][j] ,bins,minent)
				histI = create_hist(SprodI_data[k][i][j],bins,minent)
				histS = create_hist(SprodS_data[k][i][j],bins,minent)
				histS2 = create_hist(SprodS2_data[k][i][j],bins,minent)
				histA = create_hist(Action_data[k][i][j],bins,minent)
				histG = create_hist(SprodG_data[k][i][j],bins,minent)
				histG1 = create_hist(SprodG1_data[k][i][j],bins,minent)
				histG2 = create_hist(SprodG2_data[k][i][j],bins,minent)
	
				np.savetxt("/data/isilon/bause/single_particle/Sprod/hist_"+str(lag[k])+"_"+str(i)+"_"+str(j)+"_"+str(ident)+".dat", hist)
				np.savetxt("/data/isilon/bause/single_particle/Sprod/histI_"+str(lag[k])+"_"+str(i)+"_"+str(j)+"_"+str(ident)+".dat", histI)
				np.savetxt("/data/isilon/bause/single_particle/Sprod/histS_"+str(lag[k])+"_"+str(i)+"_"+str(j)+"_"+str(ident)+".dat", histS)
				np.savetxt("/data/isilon/bause/single_particle/Sprod/histS2_"+str(lag[k])+"_"+str(i)+"_"+str(j)+"_"+str(ident)+".dat", histS2)
				np.savetxt("/data/isilon/bause/single_particle/Sprod/histA_"+str(lag[k])+"_"+str(i)+"_"+str(j)+"_"+str(ident)+".dat", histA)
				np.savetxt("/data/isilon/bause/single_particle/Sprod/histG_"+str(lag[k])+"_"+str(i)+"_"+str(j)+"_"+str(ident)+".dat", histG)
				np.savetxt("/data/isilon/bause/single_particle/Sprod/histG1_"+str(lag[k])+"_"+str(i)+"_"+str(j)+"_"+str(ident)+".dat", histG1)
				np.savetxt("/data/isilon/bause/single_particle/Sprod/histG2_"+str(lag[k])+"_"+str(i)+"_"+str(j)+"_"+str(ident)+".dat", histG2)



