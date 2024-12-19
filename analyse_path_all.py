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
pos_old = pos
Sprod = 0.
Sprod2 = 0.
SprodA = 0.



sqdT = np.sqrt(dT)
sigma = np.sqrt(2*Temperature/gamma)



traj = [] 
ran_l =[]
F_l = []

checkpos = np.zeros((ms,ms))
steps = np.zeros((ms,ms))
lag_max = 101
count =0
while (count < lag_max):
	stuff = series.readline()
	if stuff.strip():
		split = stuff.split()
		xold = x
		x = float(split[0])
		ran = float(split[2])
		F= F_tanh(x,kc,U,xc) + force
		pos = get_core_pos(x,ms,shift)
		steps = steps +1
		if(pos != pos_old and pos != ms):
			pos_old = pos
			checkpos[pos,:] = 1
			steps[pos,:] = 1  
			for l in range(ms):
				if(checkpos[l,pos] ==1):
					Sprod
					steps[l,pos] =0
					checkpos[l,pos] = 0

		traj.append(x)
		ran_l.append(ran)
		F_l.append(F)		
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
SprodC_data= []
Sprod1_data= []
Sprod2_data= []
Sprod3_data= []
SprodA_data= []
SprodB_data= []
for i in range(0,ms):
    Sprod_data.append([])
    SprodC_data.append([])
    Sprod1_data.append([])
    Sprod2_data.append([])
    Sprod3_data.append([])
    SprodA_data.append([])
    SprodB_data.append([])
    for j in range(0,ms):
        Sprod_data[i].append([])
        SprodC_data[i].append([])
        Sprod1_data[i].append([])
        Sprod2_data[i].append([])
        Sprod3_data[i].append([])
        SprodA_data[i].append([])
        SprodB_data[i].append([])

#tlen = 0



steps.astype(int)

fac2 = dT/(gamma*Temperature)
fac3 = np.sqrt(8.*dT/(gamma*Temperature))
fac4 = 2.

fac1 = np.sqrt(2.*Temperature*gamma/dT) 

count =0


tracount = int(tra_len/10);

for i in range(tra_len):
	stuff = series.readline()
	if stuff.strip():
		if (i % tracount ==0):
			print(i / tracount)
		split = stuff.split()
		x = float(split[0])
		ran = float(split[2])
		F= F_tanh(x,kc,U,xc) + force
		pos = get_core_pos(x,ms,shift)
		steps = steps+1 #for all possible transitions		
		traj.append(x)
		ran_l.append(ran)
		F_l.append(F)
		pos_l.append(pos)
		if(pos != pos_old and pos != ms): 
			pos_old = pos
			checkpos[pos,:] = 1 
			steps[pos,:] =1
			checkpos[pos,pos] = 0
			for l in range(ms):
				if (checkpos[l,pos] == 1):
					Sprod1 = 0 
					Sprod2 = 0
					Sprod3 = 0
					SprodA = 0
					SprodC = 0
					length = int(steps[l,pos]) # Is there a way to force int right away?
					if (length < lag_max): 
						for j in range(1,length): # length-1 steps
							p = lag_max-length+j
							Sprod1 += fac2* F_l[p] * F_l[p]
							Sprod2 += fac3* F_l[p] * ran_l[p]
							Sprod3 += fac4* ran_l[p] * ran_l[p]
							SprodC += (traj[p] - traj[p-1])*(traj[p]-traj[p-1])/(dT*2.) + F_l[p]*dT
							SprodA += (traj[p] - traj[p-1])*(F_l[p-1]+ fac1*ran_l[p]  )/(Temperature)
							#print(F_l[p-1], fac3*ran_l[p])
						count = count+1
						#np.savetxt("/data/pckr194/bause/traj/tr1_"+str(count)+".dat",traj[lag_max-length+1:])
						Sprodg = Sprod1 + Sprod2 + Sprod3
						Sprod_data[l][pos].append(Sprodg)
						Sprod1_data[l][pos].append(Sprod1)
						Sprod2_data[l][pos].append(Sprod2)
						Sprod3_data[l][pos].append(Sprod3)
						SprodA_data[l][pos].append(SprodA/steps[l,pos])   # without #steps ? 
						SprodB_data[l][pos].append(SprodA)   
						SprodC_data[l][pos].append(SprodC)   

					
	
					steps[l,pos]=0 
					checkpos[l,pos] = 0

			
		del traj[0]
		del ran_l[0]
		del F_l[0]
		del pos_l[0]



#print("min: ", pos_min , Sprod_min, Sprod_th)
#print("max: ", pos_max , Sprod_max, Sprod_th)
#
#

bins = 100
minent = 1./(bins*3)


for i in range(ms):
	for j in range(ms):
		if Sprod_data[i][j]:
			hist  = create_hist(Sprod_data[i][j] ,bins,minent)
			hist1 = create_hist(Sprod1_data[i][j],bins,minent)
			hist2 = create_hist(Sprod2_data[i][j],bins,minent)
			hist3 = create_hist(Sprod3_data[i][j],bins,minent)
			histA = create_hist(SprodA_data[i][j],bins,minent)
			histB = create_hist(SprodB_data[i][j],bins,minent)
			histC = create_hist(SprodC_data[i][j],bins,minent)

			np.savetxt("/data/isilon/bause/single_particle/Sprod/hist_"+str(i)+"_"+str(j)+"_"+str(ident)+".dat", hist)
			np.savetxt("/data/isilon/bause/single_particle/Sprod/hist1_"+str(i)+"_"+str(j)+"_"+str(ident)+".dat", hist1)
			np.savetxt("/data/isilon/bause/single_particle/Sprod/hist2_"+str(i)+"_"+str(j)+"_"+str(ident)+".dat", hist2)
			np.savetxt("/data/isilon/bause/single_particle/Sprod/hist3_"+str(i)+"_"+str(j)+"_"+str(ident)+".dat", hist3)
			np.savetxt("/data/isilon/bause/single_particle/Sprod/histA_"+str(i)+"_"+str(j)+"_"+str(ident)+".dat", histA)
			np.savetxt("/data/isilon/bause/single_particle/Sprod/histB_"+str(i)+"_"+str(j)+"_"+str(ident)+".dat", histB)
			np.savetxt("/data/isilon/bause/single_particle/Sprod/histC_"+str(i)+"_"+str(j)+"_"+str(ident)+".dat", histC)



