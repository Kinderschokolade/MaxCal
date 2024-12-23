import numpy as np
import os
import math
import argparse
import struct

# read file into an array of binary formatted strings.
def read_binary(path):
	f = open(path,'rb')
	binlist = []
	while True:
		doub = struct.unpack('d',f.read(8))[0] 
		if not doub:
			break
	binlist.append(doub)
	
	return binlist

def read_Slist(path,ms,flag,read_pos):
	f = open(path,'rb')
	f.seek(read_pos)
	Sprod_data= []
	for i in range(0,ms):
		Sprod_data.append([])
		for j in range(0,ms):
			Sprod_data[i].append([])

	count =0
	maxbytes =  1e8
	while True:
		count = count +1
		if (count > maxbytes):
			flag = True
			read_pos = f.tell()
			break
		eightbyte = f.read(8)
		if not eightbyte:
			flag = False
			break
		n = int(struct.unpack('d',eightbyte)[0])
		m = int(struct.unpack('d',f.read(8))[0]) 
		Sprod = struct.unpack('d',f.read(8))[0] 
		Sprod_data[n][m].append(Sprod)
	return flag,read_pos,Sprod_data

def create_hist(data,bins,minent,hist_old,ran):
	hist = np.zeros((bins,2))
	if not ran: # only calc edgef, if it was not done before. Otherwise use existing ones.
		hist[:,1],edgesf = np.histogram(np.asarray(data),bins = bins, density = False)
		hist[:,1] = hist[:,1] / max( hist[:,1] )
		pos = np.squeeze(np.where(hist[:,1]>minent))
		if (np.size(pos) >1 ):
			ran = (edgesf[pos[0]], edgesf[pos[-1]+1])
	else:
		pos = np.zeros(2) #so it passes the next if statement
		
	
	if(pos.size > 1):
		hist[:,1], edgesf = np.histogram(np.asarray(data),bins = bins,range = ran, density=False )
		hist[:,1] = hist[:,1] + hist_old[:,1] # include old information
		for k in range(bins):
			hist[k,0] = (edgesf[k] + edgesf[k+1]) / 2.
	
	return ran,hist


parser = argparse.ArgumentParser()
parser.add_argument('-o', help='ident')
parser.add_argument('-l','--lag', nargs='+', help='lagtimes', required=True)
args = parser.parse_args();
ident = args.o
lagv = args.lag

lagc = len(lagv)

bins = 100;
path0 = "/data/isilon/bause/"
potential = np.loadtxt(path0+"single_particle/potential_ms_"+str(ident)+".dat")
for line in open(path0+"single_particle/param_"+str(ident)+".dat","r"):
	cont = line.split()
	if(cont[0] == 'microstates'):
		ms = int(cont[1])
	if(cont[0] == 'T'):
		Temperature = float(cont[1])
	if(cont[0] == 'extf'):
		force = float(cont[1])
	
Sp_th = np.zeros((ms,ms))	

	
ran=[]
for n in range(ms):
	ran.append([])
	for m in range(ms):
		ran[n].append([])
		if (np.abs(n-m) <15):
			Sp_th[n,m] = ( -force / ms * (n-m) + potential[n,1] - potential[m,1] ) / Temperature
		elif (m > n):
			Sp_th[n,m] = ( -force / ms * (ms+n-m)+ potential[n,1] - potential[m,1] ) / Temperature
		else:
			Sp_th[n,m] = ( +force / ms * (ms-n+m)+ potential[n,1] - potential[m,1] ) / Temperature

link =  "/usr/data/bause/single_particle/"


	
for lag in lagv:
	#Sfile = open(link+str(lag)+"/Sprod/Sprod_"+str(ident)+".dat", 'w')
	hist= np.zeros((ms,ms,bins,2))
	minent = 0.005
	flag = True
	count =0
	read_pos =0
	path = link+"Sprod_"+str(ident)+"_"+str(lag)+".dat"
	while flag:
		count = count +1
		print(count, hist[17,18,0,0], hist[17,18,-1,0])
		flag,read_pos,Slist = read_Slist(path,ms,flag,read_pos)
		for i in range(ms):
			for j in range(ms):
				if len(Slist[i][j]) > 100:
					ran[i][j],hist[i,j]  = create_hist(Slist[i][j] ,bins,minent,hist[i][j],ran[i][j])


	Smat = np.zeros((ms,ms))
	for i in range(ms):
		for j in range(ms):
			norm =  max( hist[i][j][:,1])
			if(norm > 100):
				hist[i][j][:,1] = hist[i][j][:,1] / norm
				np.savetxt("/usr/data/bause/single_particle/"+str(lag)+"/hist_"+str(i)+"_"+str(j)+"_"+str(ident)+".dat", hist[i,j])

				norm = sum(hist[i][j][:,1])
				mean = sum(hist[i][j][:,0] * hist[i][j][:,1])
				#Sfile.write(str(i) +'\t'+str(j)+'\t'+str(mean/norm)+'\t'+str(Sp_th[i,j])+'\n')
				Smat[i,j] = mean/norm
			else:
				Smat[i,j] = np.nan

	np.savetxt(link+str(lag)+"/Sprod_"+str(ident)+".dat",Smat)

	#os.system("rm "+link+"Sprod_"+str(ident)+"_"+str(lag)+".dat") # delete file at the end, because it is large 		


