import numpy as np
import os
import math
import argparse

def Gauss(x, sigma):
	return np.exp(-x[0]*x[0]/(2*sigma[0]*sigma[0])-x[1]*x[1]/(2*sigma[1]*sigma[1]));


def F_gauss(x, sigma, U, x0, lvl) :
	f = np.zeros(2)
	xshift= np.zeros(2)
	for i in range(lvl):
		for d in range(2):
			xshift[0] = x[0]-x0[i];
			xshift[1] = x[1]; # only for this case
			f[d]+= U[i] * -xshift[d]/(sigma[d]*sigma[d]) * Gauss(xshift,sigma);
	return f


def get_U_Gauss(x,sigma,U,x0,lvl):
	E = 0
	xsh = np.zeros(2)
	for j in range(lvl):
		xsh[0] = x[0] - x0[j]
		xsh[1] = x[1] # only for this case
		E -= U[j] * Gauss(xsh,sigma);
	return E 


def get_core_pos(x,ms,shift,box,MS):
	pos = [0,0]
	out =0
	for i in range(2):
		pos[i] = math.floor((x[i]/box[i]+shift[i])*ms[i]*2)
		if(pos[i]%2 == 0):
			pos[i] = pos[i]//2
		else:
			out = MS
	if (out ==0):
		out = pos[0] * ms[1] + pos[1]
	return out


parser = argparse.ArgumentParser()
parser.add_argument('-o', help='ident')
parser.add_argument('-ms1', help='end', default=30)
parser.add_argument('-ms2', help='end', default=30)
args = parser.parse_args();
ident = args.o
ms = [0,0]
ms[0] = int(args.ms1)
ms[1] = int(args.ms2)
MS = int(ms[0] * ms[1])

path = "/data/pckr194/bause/single_particle/trajectory_" +str(ident) +"_0.dat"
series = open(path, "r")
tra_len = int(os.popen("wc -l "+path).readline().split()[0])


lvl = 3
U = np.zeros(3)
xc = np.zeros(3)
xc[0] = 0.5
xc[1] = 1.5
xc[2] = 1.5

force = np.zeros(2)
sigma = np.zeros(2)

for line in open("/data/isilon/bause/single_particle/param_"+str(ident)+".dat","r"):
	cont = line.split()
	if(cont[0] == 'dT'):
		dT = float(cont[1])
	if(cont[0] == 'gamma'):
		gamma = float(cont[1])
	if(cont[0] == 'T'):
		Temperature = float(cont[1])
	if(cont[0] == 'extf0'):
		force[0] = float(cont[1])
	if(cont[0] == 'extf1'):
		force[1] = float(cont[1])
	if(cont[0] == 'U0'):
		U[0] = float(cont[1])
	if(cont[0] == 'U1'):
		U[1] = float(cont[1])
	if(cont[0] == 'U2'):
		U[2] = float(cont[1])
	if(cont[0] == 'sigma1'):
		sigma[0] = float(cont[1])
	if(cont[0] == 'sigma2'):
		sigma[1] = float(cont[1])

shift = np.zeros(2)
box = np.zeros(2)
box[0] = 3.
box[1] = 1.


x = np.zeros(2)
ran = np.zeros(2)
xold = np.zeros(2)
shift[0] = -box[0]/(ms[0]*4) # 50 % core position 
shift[1] = -box[1]/(ms[0]*4) # 50 % core position 
stuff =series.readline()
split = stuff.split()
x[0] = float(split[1])
x[1] = float(split[2])
ran[0] = float(split[3])
ran[1] = float(split[4])
pos = get_core_pos(x,ms,shift,box,MS)
pos_old = pos
Sprod = 0.
Sprod2 = 0.
SprodA = 0.

sqdT = np.sqrt(dT)


traj = [] 
ran_l =[]
F_l = []

checkpos = np.zeros((MS,MS))
steps = np.zeros((MS,MS))

lag_max = 210
lag_min = 190
count =0
while (count < lag_max):
	stuff = series.readline()
	if stuff.strip():
		split = stuff.split()
		x = np.zeros(2)
		ran = np.zeros(2)
		x = [float(split[1]), float(split[2])]
		ran = [float(split[3]), float(split[4])]
		pos = get_core_pos(x,ms,shift,box,MS)
		F= F_gauss(x,sigma,U,xc,lvl) + force
		steps = steps+1
		if (pos != pos_old and pos != MS):
			pos_old = pos
			checkpos[pos,:] =1
			steps[pos,:] =1
			for l in range(MS):
				if checkpos[l,pos]==1:
					steps[l,pos] = 0
					checkpos[l,pos] =0
		traj.append(x)
		ran_l.append(ran)
		F_l.append(F)		
		count = count +1

	
path1 = "/data/pckr194/bause/traj/tr_"

Sprod_th = np.zeros((MS,MS))
for i in range(MS):
	for j in range(MS):
            val = (i - j) / box[0] 
            if (val > box[0]/2) : 
                dist = val-box[0]
            elif (val < -box[0]/2):
                dist = box[0] + val
            elif (val == box[0]/2) or (val == -box[0]/2):
                dist = 0
            else:
                dist = val
            Sprod_th[i][j] = (-force[0]  * dist + get_U_Gauss(x,sigma,U,xc,lvl) - get_U_Gauss(x,sigma,U,xc,lvl))/Temperature #0.945 is a manual factor


Sprod_data= []
for i in range(0,MS):
    Sprod_data.append([])
    for j in range(0,MS):
        Sprod_data[i].append([])

steps.astype(int)

tra_len = 1000
count =0

si = np.sqrt(2*Temperature/gamma)
fac1 = sqdT/(gamma*si)
fac2 = dT/(2*gamma*gamma*si*si)


for i in range(tra_len):
	stuff = series.readline()
	if stuff.strip():	
		split = stuff.split()
		x = np.zeros(2)
		ran= np.zeros(2)
		x = [float(split[1]), float(split[2])]
		ran = [float(split[3]), float(split[4])]
		pos = get_core_pos(x,ms,shift,box,MS)
		F= F_gauss(x,sigma,U,xc,lvl) + force
		steps = steps+1	
		traj.append(x)
		ran_l.append(ran)
		F_l.append(F)
		# make this more effective: Most transitions never occure but are tracked anyway
		if(pos != pos_old and pos != MS):
			pos_old = pos
			checkpos[pos,:] =1
			steps[pos,:] =1
			checkpos[pos,pos] =0
			for l in range(MS):
				if (checkpos[l,pos] == 1):
					Sprod = 0 
					Sprod2 = 0
					length = int(steps[l,pos]) # Is there a way to force int right away?
					if (length < lag_max and length > lag_min):
						for j in range(1,length):
							F= F_l[lag_max-length+j]+ F_l[lag_max-j+1]
							Sprod+= fac1*np.sum(-F *ran_l[lag_max-j])
							Sprod2+= fac2*np.sum(F*F)
						Sprodg = Sprod+Sprod2
						Sprod_data[l][pos].append(Sprodg)
						count= count +1
						print(count, l, pos , Sprodg,length)
						np.savetxt(path1+str(count)+".dat",traj[lag_max-length+1:])
						steps[l,pos] = 0
						checkpos[l,pos] = 0

					if (length > lag_max):	
						steps[l,pos] = 0
						checkpos[l,pos] = 0
			
		del traj[0]
		del ran_l[0]
		del F_l[0]
		


#print("min: ", pos_min , Sprod_min, Sprod_th)
#print("max: ", pos_max , Sprod_max, Sprod_th)
#
#
bins = 100
hist= np.zeros((bins,2))
#for i in range(MS):
#	for j in range(MS):
#		if Sprod_data[i][j]:
#			
#			ra = 4* abs(Sprod_th[i,j])
#			hist[:,1],edgesf = np.histogram(np.asarray(Sprod_data[i][j]),bins = bins,range=(-ra,ra))
#			for k in range(bins):
#				hist[k,0] = (edgesf[k] + edgesf[k+1]) / 2.
#			np.savetxt("/data/isilon/bause/single_particle/Sprod/hist_"+str(i)+"_"+str(j)+"_"+str(ident)+".dat", hist)



